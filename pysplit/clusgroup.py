import numpy as np
import math
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as clrs


class ClusterGroup(object):
    """
    Class for processing and plotting Clusters
    """

    def __init__(self, cluster_object_list):
        """
        Initialize ClusterGroup object.

        Parameters
        ----------
        cluster_object_list : list of cluster objects
            Cluster objects from the same HYSPLIT clustering run.
        """

        self.totalclusters = len(cluster_object_list)
        self.clusters = cluster_object_list
        totaltraj = []
        for cluster in self.clusters:
            totaltraj.append(cluster.trajcount)
        self.totaltrajcount = sum(totaltraj)

    def set_coordinates(self, endpoints_dir):
        """
        Sets the mean coordinates in constituent clusters

        Parameters
        ----------
        endpoints_dir : string
            Full or relative path to where the c mean paths are kept.

        """

        for clus in self.clusters:

            endpoints_fname = ('C' + str(clus.clusternumber) + '_' +
                               str(self.totalclusters) + 'mean.tdump')
            endpoints_file = os.path.join(endpoints_dir, endpoints_fname)

            clus.set_coordinates(endpoints_file)

    def map_clusters(self, basemap, ax=None, figsize=(20, 20),
                     color_var='mean_mf', color_min=None, color_max=None,
                     color_rescale=None, width_var='count',
                     width_rescale=None, width_adjust=1.0,
                     colormap='blues'):
        """
        Plots mean cluster paths with color and width scaled to some variable.

        Creates a map where each cluster trajectory is plotted with a color
            determined by the mean variable value and a width by the trajectory
            count.  Color and width scaled according to user preferences.

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
        color_var : string
            Default 'mean_mf'.  ['total_raint0'|'total_mft1'|'mean_q'|
            'mean_w'|'mean_mf'|'mean_mft1'|'mean_raint0'|'total_mf']
        color_min : int or float
            Default None.  The minimum value for color mapping.  If None,
            color_min will be the minimum value of the data.
        color_max : int or float
            Default None.  The maximum value for color mapping.  If None,
            color_max will be the maximum value of the data.
        color_rescale : string
            Default None.  ['sqrt'|ln'|'log']
            Determines how color data is transformed, if at all
        width_var : string
            Default 'count'.  ['absolute_count'|'relative_count'|
            'total_raint0'|'total_mft1'|'mean_q'|
            'mean_w'|'mean_mf'|'mean_mft1'|'mean_raint0'|'total_mf']
        width_rescale : string
            Default None.  [sqrt'|'ln'|'log']
            Determines how width data is transformed, if at all
        width_adjust : float
            Default 1.0.  Used to adjust linewidths to a reasonable value
            while maintaining relationships between linesizes.
        colormap : string
            Default 'blues'.  ['jet'|'blues'|'anomaly'|'heat'|'earth']
            Passed to a dictionary which retrieves the indicated colormap

        Returns
        -------
        fig : matplotlib Figure instance, optional
            Newly created Figure instance.  Only returned if a new figure
            and axis had to be created, i.e. ax=None and `basemap`
            is a MapDesign instance.
        ax : matplotlib Axes instance, optional
            Newly created Axes instance.  Only returned if a new figure
            and axis had to be created, i.e. ax=None and `basemap`
            is a MapDesign instance.
        cavemap : Basemap instance
            Basemap instance with Trajectory moisture uptake plotted on it.
        colors : list of RGBA tuples
            The color of each cluster
        plot_order : list of ints
            The order of `colors` relative to self.clusters

        """

        try:
            if ax is None:
                fig, ax, cavemap = basemap.make_basemap(figsize)
            else:
                cavemap = basemap.make_basemap(figsize, ax=ax)
        except AttributeError:
            cavemap = basemap

        colorvar_list = []
        widthvar_list = []

        # Get color variables
        for cluster in self.clusters:
            if not hasattr(cluster, color_var):
                cluster.set_meanvar()
            colorvar_list.append(getattr(cluster, color_var))

        # Get width variables
        if 'count' in width_var:
            for cluster in self.clusters:
                widthvar_list.append(cluster.trajcount)
            if 'relative' in width_var:
                tmp = []
                for w in widthvar_list:
                    tmp.append((w / self.totaltrajcount) * 100)
                widthvar_list = tmp
        else:
            for cluster in self.clusters:
                if not hasattr(cluster, width_var):
                    cluster.set_meanvar()
                widthvar_list.append(getattr(cluster, width_var))

        if color_rescale is not None:
            c_transf = get_transform(color_rescale)

            for i in range(0, self.totalclusters):
                colorvar_list[i] = c_transf[0](colorvar_list[i])

        if width_rescale is not None:
            w_transf = get_transform(width_rescale)

            for i in range(0, self.totalclusters):
                widthvar_list[i] = ((w_transf[0](widthvar_list[i]))
                                    * width_adjust)

        if color_min is None:
            color_min = min(colorvar_list)
        # Get color max
        if color_max is None:
            color_max = max(colorvar_list)

        # Set up mapping of scalar data to RGBA values from given colormap
        cnorm = clrs.Normalize(vmin=color_min, vmax=color_max)
        scalarmap = cm.ScalarMappable(norm=cnorm, cmap=colormap)

        # Obtain index array of linewidths so thick lines plot last
        plot_order = np.argsort(widthvar_list, kind='mergesort')

        colors = []

        for i in plot_order:
            # Map values to RGBA values from given cmap, use resulting color
            color = scalarmap.to_rgba(colorvar_list[i])
            colors.append(color)
            cavemap.plot(self.clusters[i].longitude, self.clusters[i].latitude,
                         linewidth=widthvar_list[i], color=color, zorder=19,
                         latlon=True)

        try:
            return fig, ax, cavemap, colors, plot_order
        except:
            return cavemap, colors, plot_order

    def trajplot_clusterplot(self, basemap, figsize=(20, 20),
                             colors=None, orientation='horizontal',
                             cluster_lw='relative', clus_zorder=19,
                             traj_lw=3, traj_zorder=19):

        """
        Plots trajectories on one map and cluster paths on another.

        Trajectories belonging to a cluster will be plotted as lines of
            the same color, and in the cluster plot the cluster mean path
            will also be the same color.  The linewidth of the cluster mean
            path will reflect the number of member trajectories unless
            specified otherwise.

        Parameters
        ----------
        basemap : MapDesign object
            Object containing the parameters to intialize a basemap

        Keyword Arguments
        -----------------
        figsize : tuple of ints
            Default (20,20).  Dimensions of figure
        colors : list of strings or RGB tuples
            Default None.  Length of list should match self.totalclusters.
            If None, a list of random colors will be generated.
        orientation : string
            Default 'horizontal'.  ['horizontal'|'vertical'].
            Indicates whether plots are in a row or column.
        cluster_lw : int or string
            Default 'relative'.  [int|'relative'|'absolute'].
            The linewidth of the cluster mean path is set with a constant
            value, or is set with a value that depends on the relative or
            absolute number of member trajectories in the cluster
        clus_zorder : int
            Default 19.  Zorder of cluster paths
        traj_lw : int
            Default 3.  The linewdith of each trajectory.
        traj_zorder : int
            Default 19.  The zorder of trajectory paths.

        Returns
        -------
        fig : matplotlib Figure instance
            Figure with two axes
        ax_t : matplotlib Axes instance
            The axis containing the plot of trajectory paths
        ax_c : matplotlib Axes instance
            The axis containing the plot of cluster paths
        trajmap : Basemap instance
            Left/top subplot: map with trajectory paths plotted on it.
            Trajectory colors correspond to cluster
        clusmap : Basemap instance
            Right/bottom subplot: map with cluster mean paths plotted on it.
            Cluster color corresponds to that of member trajectories.
            Linewidth is given value or is the number of member trajectories.
        colors : list of RGBA tuples
            The color used for each cluster.

        """

        # Set number of rows and columns in figure
        if orientation is 'horizontal':
            row = 1
            col = 2
        else:
            row = 2
            col = 1

        # Set colors if none provided
        if colors is None:
            colors = np.random.rand(self.totalclusters, 3)
            colors = np.vsplit(colors, self.totalclusters)
            color_tmp = []
            for c in colors:
                color_tmp.append(c[0])
            colors = color_tmp

        # Set cluster mean path width by member trajectory count if not given
        if cluster_lw is 'relative':
            print 'relative'
            cluster_lw = []
            for cluster in self.clusters:
                cluster_lw.append((cluster.trajcount /
                                   float(self.totaltrajcount)) * 100.0)

        elif cluster_lw is 'absolute':
            print 'absolute'
            cluster_lw = []
            for cluster in self.clusters:
                cluster_lw.append(cluster.trajcount)

        else:
            cluster_lw = [cluster_lw]
            cluster_lw = cluster_lw * self.totalclusters

        # Make figure
        fig, (ax_t, ax_c) = plt.subplots(row, col, figsize=figsize)

        # Make a map on each axis object using self.map_prefs
        ax_t, trajmap = basemap.make_basemap(figsize, ax=ax_t)
        ax_c, clusmap = basemap.make_basemap(figsize, ax=ax_c)

        # Plot
        for cluster, color, lw in zip(self.clusters, colors, cluster_lw):
            clusmap.plot(cluster.longitude, cluster.latitude, color=color,
                         linewidth=lw, latlon=True, zorder=clus_zorder,
                         ax=ax_c)

            for traj in cluster.trajectories:
                trajmap.plot(traj.longitude, traj.latitude, color=color,
                             linewidth=traj_lw, latlon=True,
                             zorder=traj_zorder, ax=ax_t)

        return fig, ax_t, ax_c, trajmap, clusmap, colors


def scatterprep(trajgroup, variable, transform):
    """
    Gets trajectory attributes in format for plotting

    Masks and transforms specified attribute

    Parameters
    ----------
    trajgroup : TrajectoryGroup
        A TrajectoryGroup object containing trajectories with the
        specified variable
    variable : string
        The trajectory attribute to plot
    transform : string
        The transform to be applied to the data

    Keyword Arguments
    -----------------

    Returns
    -------
    data : (M) masked ndarray of floats
        The trajectory attribute, assembled in one 1D array from all
        trajectories.  Data is transformed and invalid data masked.
    lats : (M) ndarray of floats
        The latitudes of all items in data
    lons : (M) ndarray of floats
        The longitudes of all items in data

    """

    datarray = None
    latarray = None
    lonarray = None

    for traj in trajgroup.trajectories:
        dat = getattr(traj, variable)
        lat = traj.latitude
        lon = traj.longitude
        if datarray is None:
            datarray = dat
            latarray = lat
            lonarray = lon
        else:
            datarray = np.concatenate((datarray, dat))
            latarray = np.concatenate((latarray, lat))
            lonarray = np.concatenate((lonarray, lon))

    data = np.ma.masked_less_equal(datarray, -999.0)
    lons = lonarray
    lats = latarray

    if transform is not None:
        transforms = get_transform(transform)
        data = transforms[1](data)

    return data, lons, lats


def get_transform(transform):
    """
    Dictionary of transforms.  Keys are strings representing transforms,
        values are functions

    Parameters
    ----------
    transform : string
        How the data is to be transformed.  ['sqrt'|'log'|'ln']

    Returns
    -------
    transform_dict : list of functions
        The functions necessary to adjust data to specified transform

    """

    # 0 is used for cluster mapping, 1 for regular trajectory mapping
    transform_dict = {'sqrt'  : [math.sqrt, np.sqrt],
                      'log'   : [math.log10, np.log10],
                      'ln'    : [math.log, np.log]}

    transforms = transform_dict[transform]

    return transforms
