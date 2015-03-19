import numpy as np
import os
import matplotlib.pyplot as plt
import mapmaker as mm
import hyfile_handler as hh
import traj_accessory as ta
from clusgroup import ClusterGroup


class TrajectoryGroup(object):
    """
    Class for processing and plotting multiple trajectories.

    """

    def __init__(self, traj_object_list):
        """
        Initialize TrajectoryGroup object.

        Parameters
        ----------
        traj_object_list : list of trajectory objects
            Trajectories that belong in the group.

        """
        self.trajectories = traj_object_list
        self.trajcount = len(traj_object_list)
        self.directory, _ = os.path.split(traj_object_list[0].fullpath)

    def __add__(self, other):
        """
        Combine two trajectory groups into one.  Checks for duplicate
            trajectories.

        Parameters
        ----------
        other : TrajectoryGroup or Cluster
            A different TrajectoryGroup or Cluster that may or may not
            contain some of the same trajectories

        Returns
        -------
        new_self : TrajectoryGroup
            A new TrajectoryGroup from the combination of self and other

        """

        filename_ls = []
        traj_ls = []

        for traj in self.trajectories:
            filename_ls.append(traj.filename)
            traj_ls.append(traj)

        for traj in other.trajectories:
            if traj.filename not in filename_ls:
                traj_ls.append(traj)

        new_tg = TrajectoryGroup(traj_ls)

        return new_tg

    def hystats(self, variable, sort_bytime='month', iterable=True):
        """
        Gathers t=0* data so you can make your own plots.

        Parameters
        ----------
        variable : string
            Trajectory attribute that may or may not be an iterable
            (if it is not, set iterable=False)

        Keyword Arguments
        -----------------
        sort_bytime
            ['month'|'hour'|'day'|'year'|'season'|'none']
            Lets user split array of all t=0 data into a list of subarrays,
            where each subarrays represents a unique `sort_bytime`.  For more
            sophisticated sorting, or sorting by something other than
            one time division, it is recommended that the user choose
            'none' and perform the sorting themselves.
        iterable : Boolean
            Default True.  Set to False if inspecting a trajectory attribute
            that is not an iterable, like rainstatus

        Returns
        -------
        sortedvar_list : list of 1D ndarrays
            List of arrays where each array contains the var at t=0 of each
            trajectory in a  unique `sort_bytime`.
            Returned if `sort_bytime` is not 'none'.
        unique_times : list of ints
            The unique `sort_bytime`s.  Returned if `sort_bytime` is not 'none'
        var : 1D ndarray
            The array of the requested variable at t=0* for each trajectory.

        Notes
        -----
        * If var is `distance`, t=1.  If var is `total_distance`, t=-1.

        """

        # Initialize
        var = None
        time = None

        # Get the variable at each t=0
        for traj in self.trajectories:
            v = getattr(traj, variable)

            if iterable:
                if variable is 'distance':
                    v = v[1]
                elif variable is 'total_distance':
                    v = v[-1]
                else:
                    v = v[0]

            if var is None:
                var = v
            else:
                var = np.concatenate((var, v))

        # Sort data, if a `sort_bytime` is specified
        if sort_bytime is not 'none':
            for traj in self.trajectories:
                t = getattr(traj, sort_bytime)

                if sort_bytime is not 'season':
                    t = int(t[0])

                if time is None:
                    time = t
                else:
                    time = np.concatenate((time, t))

            # Find unique times
            unique_times = np.unique(var)

            # Get indices to sort by
            sorted_time = np.argsort(time, kind='mergesort')

            # Sort
            sorted_var = var[sorted_time]

            first_occurrence = []

            # Split into a list of arrays, where each array is all the data
            # in the specified time
            for u in unique_times:
                first_occurrence.append(np.nonzero(sorted_time == u))[0][0]

            sortedvar_list = np.split(sorted_var, first_occurrence[1:])

            unique_times = unique_times.tolist()
            unique_times = [int(i) for i in unique_times]

            return sortedvar_list, unique_times

        else:
            return var

    def stack_trajcoords(self):
        """
        Gathers the latitudes and longitudes of each member trajectory
            into two 1D ndarrays

        """

        trajlats = None
        trajlons = None

        for traj in self.trajectories:
            lats = traj.latitude
            lons = traj.longitude

            if trajlats is None:
                trajlats = lats
            else:
                trajlats = np.concatenate((trajlats, lats))

            if trajlons is None:
                trajlons = lons
            else:
                trajlons = np.concatenate((trajlons, lons))

        self.all_trajlats = trajlats
        self.all_trajlons = trajlons

    def grid_trajvar(self, variable, cell_value, grid_res=0.5,
                     use_wherebin=True, normalize=False):
        """
        Grids

        Parameters
        ----------
        variable : string
            Attribute to gather into one array
        cell_value : string
            Determines the value of each cell from the contents of the bin.
            ['median'|'mean'|'cumulative'|'max'|'min'|'range'|'stdev']

        Keyword Arguments
        -----------------
        grid_res : float
            Default 0.5.  The grid cell size in degrees
        use_wherebin : Boolean
            Default True.  Use wherebin to repopulate a new self.grid, if
            gridding has already occurred.  If grid_res != self.grid_res,
            then use_wherebin will be overridden.
        normalize : Boolean
            Default False.  If True, normalizes grid to a 0-1 scale.

        """

        var_array = None

        # Get stack of coordinates
        if not hasattr(self, 'all_trajlats'):
            self.stack_trajcoords()

        # Get variable stack
        for traj in self.trajectories:
            if hasattr(traj, variable):
                var = getattr(traj, variable)
                if var_array is None:
                    var_array = var
                else:
                    var_array = np.concatenate((var_array, var))
            else:
                raise AttributeError('Please set attribute and try again')

        self.var = np.ma.masked_less_equal(var_array, -999.0)

        # Grid the data if you haven't before or you don't want to use wherebin
        if (not hasattr(self, 'wherebin') or not use_wherebin
            or grid_res != self.grid_res):

            (self.grid, self.xi, self.yi,
                self.bins, self.wherebin) = ta.grid_data(self.all_trajlons,
                                                         self.all_trajlats,
                                                         self.var, cell_value,
                                                         grid_res)
            self.grid_res = grid_res
        # Grid the data using wherebin
        elif use_wherebin:

            # Initialize dictionary
            cell_value_dict = {'cumulative' : np.sum,
                               'mean' : np.mean,
                               'median' : np.median,
                               'max' : np.max,
                               'min' : np.min,
                               'stdev' : np.std,
                               'range' : ta.maxmin_diff}

            # Initialize grid
            newgrid = np.zeros(self.grid.shape, dtype=self.grid.dtype)

            # Put values into grid
            for r in range(newgrid.shape[0]):
                for c in range(newgrid.shape[1]):
                    val = self.var[self.wherebin[r][c]]
                    if val.size != 0:
                        binval = cell_value_dict[cell_value](val)
                    else:
                        binval = np.nan

                    newgrid[r, c] = binval

            # Mask invalid values
            newgrid = np.ma.masked_less_equal(newgrid, -999.0)

            # Set attribute
            self.grid = newgrid

        if normalize:
            gridmax = np.max(self.grid)
            gridmin = np.min(self.grid)

            self.grid = (self.grid - gridmin) / (gridmax - gridmin)

    def grid_moisturevar(self, uptake_type, scale, cell_value,
                         grid_res=0.5, normalize=False):
        """
        Grids moisture uptake data.

        Parameters
        ----------
        uptake_type : string
            ['both'|'above'|'below'|'all points']
        scale : string
            Moisture variable to gather into one array
            ['absolute dq'|'absolute dqi|'fractional']
        cell_value : string
            Determines the value of each cell from the contents of the bin.
            ['median'|'mean'|'cumulative'|'max'|'min'|'range'|'stdev']

        Keyword Arguments
        -----------------
        grid_res : float
            Default 0.5.  The grid cell size in degrees
        normalize : Boolean
            Default False.  If True, normalizes grid to a 0-1 scale.

        """

        # Get moisture parameters
        lon_col = self.trajectories[0].moisture_header.index('longitude')
        lat_col = self.trajectories[0].moisture_header.index('latitude')
        dq_col = self.trajectories[0].moisture_header.index('delta q')
        dqi_col = self.trajectories[0].moisture_header.index('delta q initial')
        f_col = self.trajectories[0].moisture_header.index('f')
        e_col = self.trajectories[0].moisture_header.index('e')
        d_col = self.trajectories[0].moisture_header.index('unknown fraction')

        # Uptake, scale,: data_col, data_rows, cmap
        opts = {'all points': {'fractional'  : [d_col, None],
                               'absolute dqi': [dqi_col, None],
                               'absolute dq' : [dqi_col, None]},
                'both'      : {'absolute dqi': [dqi_col, 'ef'],
                               'absolute dq' : [dq_col, 'ef']},
                'above'     : {'absolute dqi': [dqi_col, 'e'],
                               'absolute dq' : [dq_col, 'e'],
                               'fractional'  : [e_col, None]},
                'below'     : {'absolute dqi': [dqi_col, 'f'],
                               'absolute dq' : [dq_col, 'f'],
                               'fractional'  : [f_col, None]}}

        data_col, data_rows = opts[uptake_type][scale]

        moisture_data = None
        mlats = None
        mlons = None
        mdata = None

        for traj in self.trajectories:
            moistarray = traj.moisture_sources
            if moisture_data is None:
                moisture_data = moistarray
            else:
                moisture_data = np.concatenate((moisture_data, moistarray))

        if data_rows is None:
            mlons = moisture_data[:, lon_col]
            mlats = moisture_data[:, lat_col]
            mdata = moisture_data[:, data_col]
        else:
            if 'e' in data_rows:
                e_rows = np.nonzero(moisture_data[:, e_col] > -999.0)
                mlons = moisture_data[e_rows, lon_col]
                mlats = moisture_data[e_rows, lat_col]
                mdata = moisture_data[e_rows, data_col]

            if 'f' in data_rows:
                f_rows = np.nonzero(moisture_data[:, f_col] > -999.0)
                lons = moisture_data[f_rows, lon_col]
                lats = moisture_data[f_rows, lat_col]
                data = moisture_data[f_rows, data_col]

                if mlats is None:
                    mlons = lons
                    mlats = lats
                    mdata = data
                else:
                    mlons = np.concatenate((mlons, lons), axis=1)
                    mlats = np.concatenate((mlats, lats), axis=1)
                    mdata = np.concatenate((mdata, data), axis=1)

        mdata = np.ma.masked_less_equal(mdata, -999.0)

        (self.mgrid, self.mxi,
            self.myi, self.mbins, _) = ta.grid_data(mlons, mlats, mdata,
                                                    cell_value, grid_res)

        if normalize:
            gridmax = np.max(self.mgrid)
            gridmin = np.min(self.mgrid)

            self.mgrid = (self.mgrid - gridmin) / (gridmax - gridmin)

    def unique_dates(self):
        """
        Acquire unique parcel launch times, in order.

        """

        datestrs = []
        for traj in self.trajectories:
            datestrs.append(traj.datestring)

        datestrs = list(set(datestrs))

        self.datestrings = datestrs

    def set_raincount(self, reset_traj_rainstatus=False,
                      rainy_criterion='rainfall', check_steps=1,
                      rh_threshold=0.8):
        """
        Finds the number of trajectories with `rainstatus` = True.

        If the trajectories do not currently have the attribute `rainstatus` or
            the user wants to reset rainstatus (`reset_traj_rainstatus`=True),
            then the rainstatus of the trajectories will be set using the
            given criteria.

        Keyword Arguments
        -----------------
        reset_traj_rainstatus : Boolean
            Default False.  Examines preset self.trajectories `rainstatus`.
            If the trajectories do not have `rainstatus` attribute, then
            `reset_traj_rainstatus` will be automatically set to True.
            When True, `rainstatus` is set using given criteria.
        rainy_criterion : string
            ['rainfall'|'relative humidity'|'specific humidity']
            'rainfall' : set self.rainstatus to True if trajectory has
                rain within the indicated timesteps
            'relative humidity' : set self.rainstatus to True if trajectory
                is above the given rh_threshold within the indicated timesteps
            'specific humidity' : set self.rainstatus to True if trajectory
                meets the rainfall criteria and has a negative change
                in specific humidity
        check_steps : integer
            The number of timesteps from the beginning to search for rain.
        rh_threshold : float
            The relative humidity above which it is considered to be raining.):

        """

        if not hasattr(self.trajectories[0], 'rainstatus'):
            reset_traj_rainstatus = True

        if reset_traj_rainstatus:
            for traj in self.trajectories:
                traj.set_rainstatus(rainy_criterion, check_steps, rh_threshold)

        raincount = 0

        for traj in self.trajectories:
            if traj.rainstatus:
                raincount += 1

        self.raincount = raincount

    def make_infile(self):
        """
        Take trajectories in trajectory group and write path to infile

        If a specific subset of the trajectory group is needed (i.e.,
            only rainy trajectories), create a new trajectory group
            that has only trajectories meeting that criteria

        """

        infile = open(os.path.join(self.directory, 'INFILE'), 'w')

        for traj in self.trajectories:
            output = str(traj.cfullpath)
            output = output.replace('\\', '/')
            infile.writelines(output + '\n')
            infile.flush()

    def spawn_clusters(self, cfile):
        """
        Acquires the distribution of trajectories from `cfile` and
            creates new cluster objects and a cluster group object based
            on that information

        Parameters
        ----------
        cfile : string
            The filename of the CLUSLIST_# file that indicates
            trajectory distribution among clusters

        Returns
        -------
        clustergroup : ClusterGroup object
            The group of clusters derived from original TrajectoryGroup.
            A ClusterGroup consists of a list of Cluster objects, which are
            specialized TrajectoryGroup objects.

        """

        traj_inds, totalclusters = hh.load_clusterfile(cfile)

        all_clusters = []
        i = 0

        while i < totalclusters:
            cluster_number = i + 1
            trajlist = [self.trajectories[j] for j in traj_inds[i]]
            clusterobj = Cluster(trajlist, cluster_number)

            all_clusters.append(clusterobj)

            i = i + 1

        clustergroup = ClusterGroup(all_clusters)

        return clustergroup

    def map_data_line(self, basemap, ax=None, figsize=(20, 20), zorder=19,
                      show_timesteps=False, show_paths=True):
        """
        Make a plot of the trajectories in the TrajectoryGroup.

        Trajectory color and linewidth is controlled by the trajectory
            attributes of the same names.  Use trajectory.set_color() and
            trajectory.set_lw() to adjust.

        Parameters
        ----------
        basemap : MapDesign or Basemap instance
            MapDesign instances contain the parameters to initialize a Basemap.
            If a MapDesign is provided, the Basemap may be initialized on a
            given axis using the kwarg `ax`.  If ax is None, then new figure,
            axis, and Basemap instances will be created.
            A previously-generated Basemap may be provided instead of a
            MapDesign.

        Keyword Arguments
        -----------------
        ax : matplotlib axes instance
            Default None.  The axis on which to draw a new Basemap, if
            `basemap` is not a Basemap instance.  If None, and `basemap` is a
            MapDesign instance, then new figure and axis instances
            will be created.
        figsize : tuple of ints
            Default (20,20).  The size of a new figure instance, if one must be
            created.
        zorder : int
            Default 19.  The zorder of the Trajectory lines on the map.
        show_timesteps : Boolean
            Default False.  If True, marker points will be plotted.
        show_paths : Boolean
            Default True.  If False, linestyle will be none.
Returns

        -------
        fig : matplotlib Figure instance, optional
            Newly created Figure instance.  Only returned if a new figure and
            axis had to be created, i.e. ax=None and `basemap` is a
            MapDesign instance.
        ax : matplotlib Axes instance, optional
            Newly created Axes instance.  ONly returned if a new figure and
            axis had to be created, i.e. ax=None and `basemap` is a MapDesign
            instance.
        cavemap : Basemap instance
            Basemap instance with Trajectory paths plotted on it.

        """

        # Initialize basemap on given or new figure, axis
        try:
            if ax is None:
                fig, ax, cavemap = basemap.make_basemap(figsize)
            else:
                cavemap = basemap.make_basemap(figsize, ax=ax)
        except AttributeError:
            cavemap = basemap

        # Style dictionaries
        mdict = {True : 'o',
                 False : None}
        lsdict = {True : '-',
                  False : ''}

        for traj in self.trajectories:
            cavemap.plot(traj.longitude, traj.latitude, color=traj.trajcolor,
                         linewidth=traj.linewidth,
                         marker=mdict[show_timesteps], latlon=True,
                         linestyle=lsdict[show_paths], markeredgecolor='none')

        try:
            return fig, ax, cavemap
        except:
            return cavemap

    def map_data_scatter(self, basemap, variable, ax=None, figsize=(20, 20),
                         zorder=19, ptsize=25, color_min=None, color_max=None,
                         colormap='blues', alpha=1.0, rescale=None,
                         sizevar=None, size_rescale=None):
        """
        Make a scatter plot of the trajectories in the TrajectoryGroup.

        Data may be scatter plotted as color change and is rescaled or not
        according to user preferences.

        Scatter plot may have second variable plotted as size.

        Parameters
        ----------
        basemap : MapDesign or Basemap instance
            MapDesign instances contain the parameters to initialize a Basemap.
            If a MapDesign is provided, the Basemap may be initialized on a
            given axis using the kwarg `ax`.  If ax is None, then new
            figure, axis, and Basemap instances will be created.
            A previously-generated Basemap may be provided instead of a
            MapDesign.
        variable : string
            Trajectory attribute name.  The variable plotted as a color change.

        Keyword Arguments
        -----------------
        ax : matplotlib axes instance
            Default None.  The axis on which to draw a new Basemap, if
            `basemap` is not a Basemap instance.  If None, and `basemap` is a
            MapDesign instance, then new figure and axis instances
            will be created.
        figsize : tuple of ints
            Default (20,20).  The size of a new figure instance, if one must be
            created.
        zorder : int
            Default 19.  The zorder of the Trajectory lines on the map.
        ptsize : int
            Default 25.  The size of the scatter points.
        color_min : int or float
            Default None.  The minimum value for color mapping.  If None,
            color_min will be the minimum value of the data.
        color_max : int or float
            Default None.  The maximum value for color mapping.  If None,
            color_max will be the maximum value of the data.
        colormap : string
            Default 'blues'.  ['jet'|'blues'|'anomaly'|'heat'|'earth']
            The matplotlib colormap the data values are mapped to.
        alpha : float
            Default 1.0.  The opacity of the scatter points
        rescale : string
            Default None.  ['sqrt'|'ln'|'log']
            Determines how data of variable is rescaled, if at all
        sizevar : string
            Default None.  The variable to plot as a change in marker size.
        size_rescale : string
            Default None.  ['sqrt'|'ln'|'log']
            Determines how data of sizevar is rescaled, if at all

        Returns
        -------
        fig : matplotlib Figure instance, optional
            Newly created Figure instance.  Only returned if a new figure
            and axis had to be created, i.e. ax=None and `basemap`
            is a MapDesign instance.
        ax : matplotlib Axes instance, optional
            Newly created Axes instance.  ONly returned if a new figure
            and axis had to be created, i.e. ax=None and `basemap`
            is a MapDesign instance.
        cavemap : Basemap instance
            Basemap instance with Trajectory data plotted on it.
        cm : matplotlib PathCollection instance
            Mappable for use in creating colorbars.  Colorbars may be created
            using make_cbar() or make_cax_cbar().

        """

        # Initialize basemap on given or new figure, axis
        try:
            if ax is None:
                fig, ax, cavemap = basemap.make_basemap(figsize)
            else:
                cavemap = basemap.make_basemap(figsize, ax=ax)
        except AttributeError:
            cavemap = basemap

        # Gather color variable into one array, prepare data
        data, lons, lats = scatterprep(self, variable, rescale)

        if color_min is None:
            color_min = data.min()

        if color_max is None:
            color_max = data.max()

        # Gather size variable into one array, prepare data
        if sizevar is not None:
            data2 = scatterprep(self, sizevar, size_rescale)
            data2 = data2 * ptsize
            ptsize = data2

        colormap = mm.get_colormap(colormap)

        cm = cavemap.scatter(lons, lats, c=data, s=ptsize, cmap=colormap,
                             vmin=color_min, vmax=color_max, latlon=True,
                             zorder=zorder, edgecolor='none', alpha=alpha)

        try:
            return fig, ax, cavemap, cm
        except:
            return cavemap, cm

    def map_moisture(self, basemap, uptake, scale, ax=None, figsize=(20, 20),
                     zorder=20, ptsize=25, color_min=None, color_max=None,
                     alpha=1.0):
        """
        Plot moisture uptakes as a scatter plot.

        Plot all points considered for moisture uptake; plot all uptakes in
            fractional or absolute values, plot uptakes below or uptakes
            above vertical criteria in fractional or absolute values

        Parameters
        ----------
        basemap : MapDesign or Basemap instance
            MapDesign instances contain the parameters to initialize a Basemap.
            If a MapDesign is provided, the Basemap may be initialized on a
            given axis using the kwarg `ax`.  If ax is None, then new
            figure, axis, and Basemap instances will be created.
            A previously-generated Basemap may be provided instead of a
            MapDesign.
        uptake : string
            ['both'|'above'|'below'|'all points']
        scale : string
            ['absolute dq'|'absolute dqi|'fractional']

        Keyword Arguments
        -----------------
        ax : matplotlib axes instance
            Default None.  The axis on which to draw a new Basemap, if
            `basemap` is not a Basemap instance.  If None, and `basemap` is a
            MapDesign instance, then new figure and axis instances
            will be created.
        figsize : tuple of ints
            Default (20,20).  The size of a new figure instance, if one must be
            created.
        zorder : int
            Default 19.  The zorder of the Trajectory lines on the map.
        ptsize : int
            Default 25.  The size of the scatter points.
        color_min : int or float
            Default None.  The minimum value for color mapping.  If None,
            color_min will be the minimum value of the data.
        color_max : int or float
            Default None.  The maximum value for color mapping.  If None,
            color_max will be the maximum value of the data.
        alpha : float
            Default 1.0.  The opacity of the scatter points

        Returns
        -------
        fig : matplotlib Figure instance, optional
            Newly created Figure instance.  Only returned if a new figure
            and axis had to be created, i.e. ax=None and `basemap`
            is a MapDesign instance.
        ax : matplotlib Axes instance, optional
            Newly created Axes instance.  ONly returned if a new figure
            and axis had to be created, i.e. ax=None and `basemap`
            is a MapDesign instance.
        cavemap : Basemap instance
            Basemap instance with Trajectory moisture uptake plotted on it.
        cm : matplotlib PathCollection instance
            Mappable for use in creating colorbars.  Colorbars may be created
            using make_cbar() or make_cax_cbar().

        """

        try:
            if ax is None:
                fig, ax, cavemap = basemap.make_basemap(figsize)
            else:
                cavemap = basemap.make_basemap(figsize, ax=ax)
        except AttributeError:
            cavemap = basemap

        # Get moisture parameters
        lon_col = self.trajectories[0].moisture_header.index('longitude')
        lat_col = self.trajectories[0].moisture_header.index('latitude')
        dq_col = self.trajectories[0].moisture_header.index('delta q')
        dqi_col = self.trajectories[0].moisture_header.index('delta q initial')
        f_col = self.trajectories[0].moisture_header.index('f')
        e_col = self.trajectories[0].moisture_header.index('e')
        d_col = self.trajectories[0].moisture_header.index('unknown fraction')

        # Uptake, scale,: data_col, data_rows, cmap
        opts = {'all points': {'fractional'  : [d_col, None, plt.cm.Blues_r],
                               'absolute dqi': [dqi_col, None, plt.cm.Blues],
                               'absolute dq' : [dqi_col, None, plt.cm.Blues]},
                'both'      : {'absolute dqi': [dqi_col, 'ef', plt.cm.Blues],
                               'absolute dq' : [dq_col, 'ef', plt.cm.Blues]},
                'above'     : {'absolute dqi': [dqi_col, 'e', plt.cm.Blues],
                               'absolute dq' : [dq_col, 'e', plt.cm.Blues],
                               'fractional'  : [e_col, None, plt.cm.Blues]},
                'below'     : {'absolute dqi': [dqi_col, 'f', plt.cm.Blues],
                               'absolute dq' : [dq_col, 'f', plt.cm.Blues],
                               'fractional'  : [f_col, None, plt.cm.Blues]}}

        data_col, data_rows, cmap = opts[uptake][scale]

        if color_min is None:
            if 'absolute' in scale:
                color_min = 10.0
                for tr in self.trajectories:
                    datamin = np.min(tr.masked_sources[:, data_col])
                    if datamin < color_min:
                        color_min = datamin
            else:
                color_min = 0.0

        if color_max is None:
            if 'absolute' in scale:
                color_max = 0.0
                for tr in self.trajectories:
                    datamax = np.max(tr.masked_sources[:, data_col])
                    if datamax > color_max:
                        color_max = datamax
            else:
                color_max = 1.0

        for tr in self.trajectories:
            if data_rows is None:
                lons = tr.masked_sources[:, lon_col]
                lats = tr.masked_sources[:, lat_col]
                data = tr.masked_sources[:, data_col]
                cm = cavemap.scatter(lons, lats, c=data, s=ptsize, cmap=cmap,
                                     vmin=color_min, vmax=color_max,
                                     latlon=True, zorder=zorder, alpha=alpha,
                                     edgecolor='none')
            else:
                if 'e' in data_rows:
                    e_rows = np.nonzero(tr.masked_sources[:, e_col] > -999.0)
                    lons = tr.masked_sources[e_rows, lon_col]
                    lats = tr.masked_sources[e_rows, lat_col]
                    data = tr.masked_sources[e_rows, data_col]
                    cm = cavemap.scatter(lons, lats, c=data, s=ptsize,
                                         cmap=cmap, vmin=color_min,
                                         vmax=color_max, latlon=True,
                                         zorder=zorder, alpha=alpha,
                                         edgecolor='none')
                if 'f' in data_rows:
                    f_rows = np.nonzero(tr.masked_sources[:, f_col] > -999.0)
                    lons = tr.masked_sources[f_rows, lon_col]
                    lats = tr.masked_sources[f_rows, lat_col]
                    data = tr.masked_sources[f_rows, data_col]
                    cm = cavemap.scatter(lons, lats, c=data, s=ptsize,
                                         cmap=cmap, vmin=color_min,
                                         vmax=color_max, latlon=True,
                                         zorder=zorder, alpha=alpha,
                                         edgecolor='none')
        try:
            return fig, ax, cavemap, cm
        except:
            return cavemap, cm

    def gridmap(self, basemap, ax=None, figsize=(20, 20), ismoisture=False,
                usecontourf=False, mapcount=False, color_min=None,
                color_max=None, colormap='blues', zorder=20):
        """
        Create a pcolor or filled contourplot of gridded data.

        Parameters
        ----------
        basemap : MapDesign or Basemap instance
            MapDesign instances contain the parameters to initialize a Basemap.
            If a MapDesign is provided, the Basemap may be initialized on a
            given axis using the kwarg `ax`.  If ax is None, then new
            figure, axis, and Basemap instances will be created.
            A previously-generated Basemap may be provided instead of a
            MapDesign.

        Keyword Arguments
        -----------------
        ax : matplotlib axes instance
            Default None.  The axis on which to draw a new Basemap, if
            `basemap` is not a Basemap instance.  If None, and `basemap` is a
            MapDesign instance, then new figure and axis instances
            will be created.
        figsize : tuple of ints
            Default (20,20).  The size of a new figure instance, if one must be
            created.
        ismoisture : Boolean
            Default False.  Access gridded moisture uptake data if True.
            If False, access regular gridded data.
        usecontourf : Boolean
            Default False.  Use contourf (True) or pcolormesh (False) to
            plot gridded data.
        mapcount : Boolean
            Default False.  If true, plot counts in each gridbox.  Else,
            plot gridded variable.
        color_min : int or float
            Default None.  The minimum value for color mapping.  If None,
            color_min will be the minimum value of the data.
        color_max : int or float
            Default None.  The maximum value for color mapping.  If None,
            color_max will be the maximum value of the data.
        colormap : string
            Default 'blues'.  ['jet'|'blues'|'anomaly'|'heat'|'earth']
            Passed to a dictionary which retrieves the indicated colormap
        zorder : int
            Default 20.  The zorder of the gridded data

        Returns
        -------
        fig : matplotlib Figure instance, optional
            Newly created Figure instance.  Only returned if a new figure
            and axis had to be created, i.e. ax=None and `basemap`
            is a MapDesign instance.
        ax : matplotlib Axes instance, optional
            Newly created Axes instance.  ONly returned if a new figure
            and axis had to be created, i.e. ax=None and `basemap`
            is a MapDesign instance.
        cavemap : Basemap instance
            Basemap instance with Trajectory moisture uptake plotted on it.
        cm : matplotlib PathCollection instance
            Mappable for use in creating colorbars.  Colorbars may be created
            using make_cbar() or make_cax_cbar().

        """

        try:
            if ax is None:
                fig, ax, cavemap = basemap.make_basemap(figsize)
            else:
                cavemap = basemap.make_basemap(figsize, ax=ax)
        except AttributeError:
            cavemap = basemap

        # Prepare grids for moisture variable or other variable
        if ismoisture:
            x = self.mxi
            y = self.myi

            if mapcount:
                data = self.mbins
                data = np.ma.masked_equal(data, 0)
            else:
                data = self.mgrid

        else:
            x = self.xi
            y = self.yi

            if mapcount:
                data = self.bins
                data = np.ma.masked_equal(data, 0)
            else:
                data = self.grid

        if color_min is None:
            color_min = np.min(data)
        if color_max is None:
            color_max = np.max(data)

        colormap = mm.get_colormap(colormap)

        if usecontourf:
            cm = cavemap.contourf(x, y, data,
                                  np.linspace(color_min, color_max, num=11),
                                  cmap=colormap, latlon=True, zorder=zorder)
        else:
            cm = cavemap.pcolormesh(x, y, data, latlon=True, cmap=colormap,
                                    vmin=color_min, vmax=color_max,
                                    zorder=zorder)

        try:
            return fig, ax, cavemap, cm
        except:
            return cavemap, cm


class Cluster(TrajectoryGroup):
    """
    A special subclass of TrajectoryGroup for trajectories that have been
        clustered together using HYSPLIT's clustering process.
        Contains TrajectoryGroup attributes and functions, but also has
        Trajectory-like attributes and functions associated with it, since
        a Cluster may be represented as a mean trajectory
    """

    def __init__(self, traj_object_list, cluster_number):
        """
        Initialize Cluster object.

        Parameters
        ----------
        traj_object_list : list of trajectory objects
            Trajectories that belong in the cluster.
        cluster_number : int
            The Cluster identification number.  Distinguishes Cluster
            from other Clusters in its ClusterGroup

        """
        TrajectoryGroup.__init__(self, traj_object_list)
        self.start_longitude = traj_object_list[0].longitude[0]
        self.clusternumber = cluster_number

    def __add__(self, other):
        """
        Prints notice before calling TrajectoryGroup.__add__()

        """

        print "Basic TrajectoryGroup created, cluster methods unavailable"

        new_tg = TrajectoryGroup.__add__(self, other)

        return new_tg

    def set_coordinates(self, endpoints_file):
        """
        Initialize the coordinates of the Cluster mean trajectory path.

        Parameters
        ----------
        endpoints_dir : string
            Full or relative path to the HYSPLIT file containing the
            trajectory path

        """

        data, header, _ = hh.load_hysplitfile(endpoints_file)

        self.latitude = data[0][:, header.index('Latitude')]
        self.longitude = data[0][:, header.index('Longitude')]

        # Longitudes should be -180 to 180
        for lon in self.longitude:
            if lon > 180.0:
                lon = lon - 360.0

    def set_vector(self):
        """
        Calculate mean bearing of Cluster path and bearings between timesteps

        """

        # Coordinates must be loaded first
        if not hasattr(self, 'latitude'):
            raise ValueError('self.latitude and self.longitude missing \n' +
                             'perform set_coordinates() first')

        self.meanvector, self.bearings = ta.tracemean_vector(self.latitude,
                                                             self.longitude)

    def set_distance(self):
        """
        Calculate the distance between timesteps fo the Cluster path and the
            cumulative distance at each time step

        """

        # Coordinates must be loaded first
        if not hasattr(self, 'latitude'):
            raise ValueError('self.latitude and self.longitude missing \n' +
                             'perform set_coordinates() first')
        self.distance = ta.distance_overearth(self.latitude, self.longitude)
        self.total_distance = ta.sum_distance(self.distance)

    def set_meanvar(self):
        """
        Calculate means and sums of several moisture-related variables.

        Means include the mean of every value along every trajectory in the
            Cluster and the mean of every trajectory's latest value.  Sums are
            the sum of every trajectory's latest value.

        Attributes set:
            mean_mf
            mean_w
            mean_q
            mean_raint0
            total_raint0
            mean_mft1
            total_mft1

        """

        mf_list = []
        q_list = []
        w_list = []
        raint0_list = []
        mft1_list = []

        for traj in self.trajectories:
            # Get trajectory attributes, if necessary setting them first
            if not hasattr(traj, 'specific_humidity'):
                traj.set_specifichumidity()

            if not hasattr(traj, 'mixing_ratio'):
                traj.set_mixingratio()

            if not hasattr(traj, 'moistureflux'):
                traj.calculate_moistureflux()

            mf_list.extend(traj.moistureflux[:-1].tolist())
            w_list.extend(traj.mixing_ratio.tolist())
            q_list.extend(traj.specific_humidity.tolist())
            raint0_list.append(traj.rainfall[0])
            mft1_list.append(traj.moistureflux[0])

        self.mean_mf = np.mean(mf_list)
        self.mean_w = np.mean(w_list)
        self.mean_q = np.mean(q_list)

        self.mean_raint0 = np.mean(raint0_list)
        self.total_raint0 = np.sum(raint0_list)
        self.mean_mft1 = np.mean(mft1_list)
        self.total_mft1 = np.sum(mft1_list)
        self.total_mf = np.sum(mf_list)
