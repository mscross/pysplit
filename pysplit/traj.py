import numpy as np
import os
import hyfile_handler as hh
import traj_accessory as ta
import mapmaker as mm


class Trajectory:
    """
    Class for processing individual HYSPLIT back trajectories.

    """

    def __init__(self, trajdata, trajheader, path):
        """
        Initialize (back) trajectory object.

        Parameters
        ----------
        trajdata : (M, N) ndarray of floats
            The data array corresponding to a single HYSPLIT
            back trajectory
        trajheader : list of N strings
            Column headers for trajdata
        path : string
            The full path to the original HYSPLIT file

        """

        season_dict = {12: 'winter', 1 : 'winter', 2 : 'winter',
                       3 : 'spring', 4 : 'spring', 5 : 'spring',
                       6 : 'summer', 7 : 'summer', 8 : 'summer',
                       9 : 'autumn', 10: 'autumn', 11: 'autumn'}

        self.data = trajdata
        self.header = trajheader
        self.fullpath = path
        self.folder, self.filename = os.path.split(self.fullpath)

        if os.path.isdir(os.path.join(self.folder, 'clippedtraj')):
            self.cfolder = os.path.join(self.folder, 'clippedtraj')
            if os.path.exists(os.path.join(self.cfolder,
                                           self.filename + 'CLIPPED')):
                self.cfilename = self.filename + 'CLIPPED'
                self.cfullpath = os.path.join(self.cfolder, self.cfilename)

        self.latitude = self.data[:, self.header.index('Latitude')]
        self.longitude = self.data[:, self.header.index('Longitude')]
        self.pressure = self.data[:, self.header.index('PRESSURE')]
        self.altitude = self.data[:, self.header.index('Altitude (magl)')]
        self.temperature = self.data[:, self.header.index('AIR_TEMP')]
        self.rainfall = self.data[:, self.header.index('RAINFALL')]
        self.mixdepth = self.data[:, self.header.index('MIXDEPTH')]

        self.year = self.data[:, self.header.index('Year')]
        self.month = self.data[:, self.header.index('Month')]
        self.day = self.data[:, self.header.index('Date')]
        self.hour = self.data[:, self.header.index('Hour (UTC)')]

        self.timesteps = self.data[:, self.header.index('Time step (hr)')]
        self.sim_length = self.data.shape[0] - 1
        self.datestring = os.path.split(path)[1][-8:]

        self.season_t0 = season_dict[self.month[0]]

        self.trajcolor = 'black'
        self.linewidth = 2

    def set_rainstatus(self, rainy_criterion='rainfall', check_steps=1,
                       rh_threshold=0.8):
        """
        Determines if a trajectory is 'rainy', sets self.rainstatus
            accordingly.

        The determination is made by examining the user-given criteria
            within the user-given window of timesteps.

        Keyword Arguments
        -----------------
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
            The relative humidity above which it is considered to be raining.

        """

        is_rainy = False

        if rainy_criterion == 'relative humidity':
            # Check that relative humidity data is available
            if not hasattr(self, 'relative_humidity'):
                self.set_relativehumidity()
            if np.any(self.relativehumidity[:check_steps]) > rh_threshold:
                is_rainy = True

        else:
            # Test if rainfall is recorded
            if np.any(self.rainfall[:check_steps] > 0.0):
                if rainy_criterion == 'rainfall':
                    is_rainy = True
                else:
                    # Further refine results by examining
                    # examining changes in specific humidity
                    if not hasattr(self, 'specific_humidity'):
                        self.set_specifichumidity()

                    q0 = self.specific_humidity[0]
                    qx = self.specific_humidity[check_steps * 2]
                    dq = q0 - qx

                    if dq < 0:
                        is_rainy = True

        self.rainstatus = is_rainy

    def set_trajcolor(self, color=None):
        """
        Set the color of the trajectory path as it will appear in
            map_data_line()

        Keyword Arguments
        -----------------
        color : string or tuple of floats
            Default None.  If None, self.trajcolor will be reset to black.

        """

        if color is None:
            self.trajcolor = 'black'
        else:
            self.trajcolor = color

    def set_linewidth(self, lw=None):
        """
        Set the linewidth of the trajectory path as it will appear in
            map_data_line()

        Keyword Arguments
        -----------------
        lw : int or float
            Default None.  If None, self.linewidth will be reset to 2.
        """

        if lw is None:
            self.linewidth = 2
        else:
            self.linewidth = lw

    def set_vector(self):
        """
        Calculate mean bearing of trajectory and bearings between timesteps

        """

        self.meanvector, self.bearings = ta.tracemean_vector(self.latitude,
                                                             self.longitude)

    def set_relativehumidity(self):
        """
        Acquire relative humidity data for each timestep
            from self.data or by calculating relative humidity from other data

        """

        if 'RELHUMID' in self.header:
            self.relative_humidity = self.data[:,
                                               self.header.index('RELHUMID')]
        else:
            if hasattr(self, 'mixing_ratio'):
                self.relative_humidity = ta.convert_w2rh(self.mixing_ratio,
                                                         self.temperature,
                                                         self.pressure)
            else:
                self.set_mixingratio()
                self.relative_humidity = ta.convert_w2rh(self.mixing_ratio,
                                                         self.temperature,
                                                         self.pressure)

    def set_mixingratio(self):
        """
        Acquire mixing ratio data for each timestep
            from self.data or by calculating mixing ratio from other data

        """

        if 'H2OMIXRA' in self.header:
            self.mixing_ratio = self.data[:, self.header.index('H2OMIXRA')]
        else:
            if hasattr(self, 'specific_humidity'):
                self.mixing_ratio = ta.convert_q2w(self.specific_humidity)
            elif 'SPCHUMID' in self.header:
                self.set_specifichumidity()
                self.mixing_ratio = ta.convert_q2w(self.specific_humidity)

            elif hasattr(self, 'relative_humidity'):
                self.mixing_ratio = ta.convert_rh2w(self.relative_humidity,
                                                    self.temperature,
                                                    self.pressure)
            elif 'RELHUMID' in self.header:
                self.set_relativehumidity()

            else:
                raise AttributeError('Not enough information to ' +
                                     'calculate mixing ratio!')

    def set_specifichumidity(self):
        """
        Acquire specific humidity data for each timestep
            from self.data or by calculating specific humidity from other data

        """

        if 'SPCHUMID' in self.header:
            self.specific_humidity = self.data[:,
                                               self.header.index('SPCHUMID')]
        else:
            if hasattr(self, 'mixing_ratio'):
                self.specific_humidity = ta.convert_w2q(self.mixing_ratio)
            elif ('H2OMIXRA' in self.header or
                  hasattr(self, 'relative_humidity')):
                self.set_mixingratio()
                self.specific_humidity = ta.convert_w2q(self.mixing_ratio)

            elif 'RELHUMID' in self.header:
                self.set_relativehumidity()
                self.set_mixingratio()
                self.ta.convert_w2q(self.mixing_ratio)

            else:
                raise AttributeError('Not enough information to ' +
                                     'calculate specific humidity!')

    def set_distance(self):
        """
        Calculate the distances between each timestep and the distance between
            each timestep and time zero.

        """

        self.distance = ta.distance_overearth(self.latitude, self.longitude)

        self.total_distance = ta.sum_distance(self.distance)

    def dq_dw_drh(self):
        """
        Calculate the change in moisture between each timestep

        Works for specific humidity, mixing ratio, and relative humidity.

        """

        if not hasattr(self, 'specific_humidity'):
            self.set_specifichumidity()
        if not hasattr(self, 'relative_humidity'):
            self.set_relativehumidity()
        if not hasattr(self, 'mixing_ratio'):
            self.set_mixingratio()

        dq_list = []
        dw_list = []
        drh_list = []

        for i in range(0, self.sim_length):
            dq = self.specific_humidity[i] - self.specific_humidity[i + 1]
            dw = self.mixing_ratio[i] - self.mixing_ratio[i + 1]
            drh = self.relative_humidity[i] = self.relative_humidity[i + 1]

            dq_list.append(dq)
            dw_list.append(dw)
            drh_list.append(drh)

        self.dq = np.asarray(dq_list).astype(np.float64)
        self.dw = np.asarray(dw_list).astype(np.float64)
        self.drh = np.asarray(drh_list).astype(np.float64)

        self.dq = np.pad(self.dq, (0, 1), 'constant',
                         constant_values=(-999.0, -999.0))
        self.dw = np.pad(self.dw, (0, 1), 'constant',
                         constant_values=(-999.0, -999.0))
        self.drh = np.pad(self.drh, (0, 1), 'constant',
                          constant_values=(-999.0, -999.0))

        self.dq = np.ma.masked_less_equal(self.dq, -999.0)
        self.dw = np.ma.masked_less_equal(self.dw, -999.0)
        self.drh = np.ma.masked_less_equal(self.drh, -999.0)

    def calculate_moistureflux(self, qtype='specific_humidity'):
        """
        Calculate the moisture flux between each timestep
            using distance and the chosen humidity parameter

        Keyword Arguments
        -----------------
        qtype : string
            ['specific_humidity'|'mixing_ratio']
            Default 'specific_humidity'.

        """

        if not hasattr(self, qtype):
            if qtype == 'specific_humidity':
                self.set_specifichumidity()
            elif qtype == 'mixing_ratio' :
                self.set_mixingratio()

        if not hasattr(self, 'Distance'):
            self.set_distance()

        if qtype == 'specific_humidity':
            moisture = self.specific_humidity
        elif qtype == 'mixing_ratio':
            moisture = self.mixing_ratio

        speed = self.distance / 3600
        mf = speed[1:] * moisture[:-1]

        mf = np.pad(mf, (0, 1), 'constant', constant_values=(-999.0, -999.0))
        self.moistureflux = np.ma.masked_less_equal(mf, -999.0)

    def last_moistureuptake(self, q_type='specific_humidity'):
        """
        Method for finding moisture sources according to Gustafsson 2010.

        There are three basic scenarios for finding moisure sources.
            In the first, the trajectory follows the case outlined in
            Gustafsson 2010, where the trajectory experiences a decrease
            in q in the target region- the source region is easily
            identifiable.  The second is where the parcel shows an increase
            in q to the target area
        """

    def moistureuptake(self, rainout_threshold, evap_threshold,
                       uptake_window=6, window_overlap=0,
                       vertical_criterion='pbl', pressure_threshold=900.0,
                       mixdepth_factor=1, q_type='specific_humidity'):

        """
        Calculate the moisture flux between each timestep
            using distance and the chosen humidity parameter

        Parameters
        ----------
        rainout_threshold : float
            The level below which changes in humidity are considered to
            indicate loss due to precipitation.  Recommended -0.2 or lower
        evap_threshold : float
            The level above which changes in humidity are considered to
            indicate moisture uptake.  Recommended 0.2 to 0.5

        Keyword Arguments
        -----------------
        uptake_window : int
            Default 6.  The length in hours of the period between q0 and qx.
            Assumes that either evaporation or precipitation dominate over
            a short period of time, such as 6 hours
        window_overlap : int
            Default 0.  The number of hours each uptake window should overlap
            the previous one.  A value of 0 indicates that uptake windows
            should be discrete.  The latest timepoint in the earliest window
            will still start at endpoint + uptake_window hours.
        vertical_criterion : string
            Default 'pbl'.  ['pbl'|'prs'|'both']
            Criterion for differentiating surficial and other sources of
            moisture.
        pressure_threhold : float
            Default 900.0. The pressure level defined as equivalent
            to the planetary boundary layer (parcel is considered in
            communication with sample )
        mixdepth_factor : int or float
            Default 1.  The value by which to adjust the mixed layer depth.
            Use if mixed layer depth seems to be under- or over-estimated
        q_type : string
            Default 'specific_humidity'.
            ['specific_humidity'|'mixing_ratio']
            The humidity parameter to inspect for moisture changes

        """

        # Initialize moisture_header and empty array
        moisture_header = ['year',
                           'month',
                           'day',
                           'hour',
                           'timesteps',
                           'latitude',
                           'longitude',
                           'total_distance',
                           'average pressure',
                           'average mixDepth',
                           'average altitude',
                           'q',
                           'delta q initial',
                           'delta q',
                           'e',
                           'f',
                           'unknown fraction',
                           'total e',
                           'total f']

        # Initialize indices from moisture_array
        f_index = moisture_header.index('f')
        e_index = moisture_header.index('e')
        dq_index = moisture_header.index('delta q')
        q_index = moisture_header.index('q')

        d_total_index = moisture_header.index('unknown fraction')
        e_total_index = moisture_header.index('total e')
        f_total_index = moisture_header.index('total f')

        # Initialize empty array
        moisture_array = np.empty((0, len(moisture_header))).astype(np.float64)

        # Build the initial timepoint
        earliest = self.timesteps.size - 1
        # print earliest
        initial_timept = []

        for item in moisture_header[:moisture_header.index('total_distance')
                                    + 1]:
            initial_timept.append(getattr(self, item)[-1])

        initial_timept.extend([self.pressure[-1], self.mixdepth[-1],
                               self.altitude[-1], getattr(self, q_type)[-1]])

        # Initialize dqi, dq, e, f, dtot, etot, ftot
        initial_timept.extend([-999, -999, -999, -999, 1.0, 0.0, 0.0])

        # Put the initial timepoint into the moisture uptake array
        initial_timept = np.asarray(initial_timept).astype(np.float64)
        initial_timept = np.atleast_2d(initial_timept)

        moisture_array = np.concatenate((initial_timept,
                                         moisture_array), axis=0)

        # Initialize a list of the last timestep of each window
        # Spacing between these are determined by the uptake_window and overlap
        # print earliest - (uptake_window+1)
        # print uptake_window - window_overlap
        timesteps = range(0, earliest - (uptake_window - 1),
                          uptake_window - window_overlap)[::-1]
        print timesteps

        # Cycle through each chunk of time
        for ts in timesteps:

            timepoint = []

            # Get year, month, day, hour, and timestep at the latest
            # point in ts
            for item in moisture_header[:moisture_header.index('timesteps')
                                        + 1]:
                timepoint.append(getattr(self, item)[ts])

            # Get the latitude, longitude, and Total Distance
            # at the midpoint between the latest point in ts and the latest
            # point of the previous timestep
            if uptake_window % 2 == 0:
                midpoint = uptake_window / 2
            elif uptake_window == 1:
                midpoint = 0
            else:
                midpoint = (uptake_window + 1) / 2

            for item in moisture_header[moisture_header.index('latitude'):
                                        moisture_header.index('total_distance')
                                        + 1]:
                timepoint.append(getattr(self, item)[ts + midpoint])

            # Get the average pressure, mixdepth, and altitude over the window
            pressure = np.mean(self.pressure[ts:ts + uptake_window])
            mixdepth = (np.mean(self.mixdepth[ts:ts + uptake_window])
                        * mixdepth_factor)
            altitude = np.mean(self.altitude[ts:ts + uptake_window])

            # Find q and initial dq from previous q
            q = getattr(self, q_type)[ts]
            dqi = q - moisture_array[0, q_index]

            timepoint.extend([pressure, mixdepth, altitude, q, dqi])

            # Initialize vertical criteria
            if vertical_criterion is 'prs':
                if pressure > pressure_threshold:
                    below_vert_criteria = True
                else:
                    below_vert_criteria = False
            elif vertical_criterion is 'pbl':
                if altitude - mixdepth < 0:
                    below_vert_criteria = True
                else:
                    below_vert_criteria = False
            elif vertical_criterion is 'both':
                if (altitude - mixdepth < 0) and (pressure >
                                                  pressure_threshold):
                    below_vert_criteria = True
                else:
                    below_vert_criteria = False

            # If moisture uptake has occurred:
            if dqi > evap_threshold:

                dq = dqi

                # Adjust previous f and e fractions
                for i in range(0, moisture_array.shape[0]):

                    if moisture_array[i, e_index] > -999.0:
                        mostrecent_dq = moisture_array[i, dq_index]
                        updated_frac = mostrecent_dq / q
                        moisture_array[i, e_index] = updated_frac

                    elif moisture_array[i, f_index] > -999.0:
                        mostrecent_dq = moisture_array[i, dq_index]
                        updated_frac = mostrecent_dq / q
                        moisture_array[i, f_index] = updated_frac

                # Initialize current values of e and f
                # Find previous values of dq where there is e and f
                # Get new e_total, f_total
                if below_vert_criteria:
                    e       = -999.0
                    e_inds  = np.nonzero(moisture_array[:, e_index] > -999.0)
                    e_total = np.sum(moisture_array[e_inds, dq_index]) / q

                    f       = dq / q
                    f_inds  = np.nonzero(moisture_array[:, f_index] > -999.0)
                    f_total = (np.sum(moisture_array[f_inds, dq_index])
                               + dq) / q

                else:
                    e       = dq / q
                    e_inds  = np.nonzero(moisture_array[:, e_index] > -999.0)
                    e_total = (np.sum(moisture_array[e_inds, dq_index])
                               + dq) / q

                    f       = -999.0
                    f_inds  = np.nonzero(moisture_array[:, f_index] > -999.0)
                    f_total = np.sum(moisture_array[f_inds, dq_index]) / q

                d_total = 1.0 - (e_total + f_total)

            # If precipitation has occurred:
            elif dqi < rainout_threshold:

                dq = -999

                # Adjust previous dq values
                for i in range(0, moisture_array.shape[0]):

                    if moisture_array[i, f_index] > -999.0:
                        mostrecent_frac = moisture_array[i, f_index]
                        updated_dq = mostrecent_frac * q
                        moisture_array[i, dq_index] = updated_dq

                    elif moisture_array[i, e_index] > -999.0:
                        mostrecent_frac = moisture_array[i, e_index]
                        updated_dq = mostrecent_frac * q
                        moisture_array[i, dq_index] = updated_dq

                # Initialize current e and f
                e = -999
                f = -999

                # Copy previous total fractions
                e_total = moisture_array[0, e_total_index]
                f_total = moisture_array[0, f_total_index]
                d_total = moisture_array[0, d_total_index]

            # If dqi falls between the rainout and evaporative thresholds:
            else:

                # Initialize
                dq = -999
                e  = -999
                f  = -999

                # Copy previous values
                e_total = moisture_array[0, e_total_index]
                f_total = moisture_array[0, f_total_index]
                d_total = moisture_array[0, d_total_index]

            # Add newest timestep to moisture_array
            timepoint.extend([dq, e, f, d_total, e_total, f_total])
            timepoint = np.atleast_2d(np.asarray(timepoint).astype(np.float64))

            moisture_array = np.concatenate((timepoint, moisture_array),
                                            axis=0)

        # Mask missing/invalid data
        masked_moistarr = np.ma.masked_less_equal(moisture_array, -999)

        self.moisture_sources = moisture_array
        self.masked_sources = masked_moistarr
        self.moisture_header = moisture_header

    def load_forwardtraj(self, forward_dir):
        """
        Acquires data from forward trajectory launched from earliest timepoint.

        Parameters
        ----------
        forward_dir : string
            Location of forward trajectory (do not include filename)

        """

        # Check for directory
        if not os.path.isdir(forward_dir):
            raise OSError('Forward trajectory directory does not exist!')

        # Construct filename of forward trajectory
        self.forwardpath = os.path.join(forward_dir, self.filename + 'FORWARD')

        # Check that this trajectory has a corresponding forward trajectory
        if not os.path.exists(self.forwardpath):
            raise OSError('File not found: ' + self.forwardpath)

        fdatlist, _, _ = hh.load_hysplitfile(self.forwardpath)

        # Check that there is only one trajectory in the file
        if len(fdatlist) > 1:
            print 'Multiple-trajectory files not supported!'
            self.forwardpath = None
            self.fdata = None
        else:
            self.fdata = fdatlist[0]
            self.flatitude = self.fdata[:, self.header.index('Latitude')]
            self.flongitude = self.fdata[:, self.header.index('Longitude')]
            self.faltitude = self.fdata[:,
                                        self.header.index('Altitude (magl)')]
            self.fdistance = ta.distance_overearth(self.flatitude,
                                                   self.flongitude)
            self.ftotal_dist = ta.sum_distance(self.fdistance)

    def integration_error(self):
        """
        Estimates the xy integration error based on distance between
            back/forward start/end points and the total distance traveled.
            Also estimates the z integration error based on altitude difference
            between back/forward start/end points and the total range of
            altitude encountered
        Integration error is a minor source of trajectory error.  Resolution
            error should also be checked.

        """

        # Check for necessary data
        if not hasattr(self, 'ftotal_dist'):
            raise AttributeError('Forward trajectory data missing!')
        if not hasattr(self, 'total_distance'):
            self.set_distance()

        # Acquire total horizontal and vertical distances
        f_distance = self.ftotal_dist[-1]
        b_distance = self.total_distance[-1]
        falt_range = np.max(self.faltitude) - np.min(self.faltitude)
        balt_range = np.max(self.altitude) - np.min(self.altitude)

        # Gather coordinates of back trajectory launch location
        # according to both back and forward trajectories.
        site_lats = [self.latitude[0], self.flatitude[-1]]
        site_lons = [self.longitude[0], self.flongitude[-1]]

        site_lats = np.asarray(site_lats).astype(np.float64)
        site_lons = np.asarray(site_lons).astype(np.float64)

        # Horizontal distance between back traj launch and forward traj end pts
        site_distance = ta.distance_overearth(site_lats, site_lons)[1]

        # Vertical distance between back traj launch and forward traj end pts
        z_distance = self.altitude[0] - self.faltitude[-1]

        # Distance between points divided by total h or v distance
        self.integ_error_xy = (((site_distance / (f_distance + b_distance))
                               * 100) / 2)
        self.integ_error_z = (((z_distance / (falt_range + balt_range))
                              * 100) / 2)

    def map_traj_scatter(self, cavemap, variable, sizevar=None, **kwargs):
        """
        Scatter plot of trajectory data

        Parameters
        ----------
        cavemap : Basemap instance
            Initialize a basemap first using MapDesign.make_basemap()
        variable : string
            The variable to plot as a color change

        Keyword Arguments
        -----------------
        sizevar : string
            Default None.  The variable to plot as a marker size change

        Other Parameters
        ----------------
        kwargs : passed to traj_scatter() and then ax.scatter()

        Returns
        -------
        cm : matplotlib PathCollection instance
            Mappable for use in creating colorbars.  Colorbars may be created
            in PySPLIT using make_cbar() or make_cax_cbar()

        """

        data = getattr(self, variable)

        if sizevar is not None:
            sizedata = getattr(self, variable)
        else:
            sizedata = None

        cm = mm.traj_scatter(data, self.longitude, self.latitude,
                             cavemap, sizedata=sizedata, **kwargs)

        return cm

    def map_traj_path(self, cavemap, color=None, lw=None, **kwargs):
        """
        Line plot of trajectory path.

        Parameters
        ----------
        cavemap : Basemap instance
            Initialize a basemap first using MapDesign.make_basemap()

        Keyword Arguments
        -----------------
        color : string, tuple
            Default None.  If None, then traj.color will be used.  Any
            matplotlib-approved color accepted
        lw : int
            Default None.  If None, then traj.linewidth will be used.

        Other Parameters
        ----------------
        kwargs passed to traj_scatter() and then Basemap.plot(),
            and then axis.plot()

        """

        if color is None:
            color = self.trajcolor
        if lw is None:
            lw = self.linewidth

        mm.traj_path(cavemap, self.longitude, self.latitude,
                     color, lw, **kwargs)
