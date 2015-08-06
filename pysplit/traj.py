from __future__ import division
import numpy as np
import pandas as pd
import os
import hyfile_handler as hh
import traj_accessory as ta
import mapmaker as mm


class Trajectory(pd.DataFrame):
    """
    Class for processing individual HYSPLIT back trajectories.
    Subclass of pandas DataFrame.

    """

    def __init__(self, trajdata, datetime, trajheader, folder, filename,
                 cfolder=None):
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
        name = ''
        if trajdata[0,0] < 1:
            name = '_' + trajdata[0,0]

        pd.DataFrame.__init__(self, data=trajdata[:, 1:], columns=trajheader[1:])

        self['DateTime'] = datetime
        self.set_index('Timestep')

        self.rename(columns={'AIR_TEMP' : 'Temperature',
                             'PRESSURE' : 'Pressure',
                             'RAINFALL' : 'Rainfall',
                             'MIXDEPTH' : 'Mixing_Depth',
                             'RELHUMID' : 'Relative_Humidity',
                             'H2OMIXRA' : 'Mixing_Ratio',
                             'SPCHUMID' : 'Specific_Humidity',
                             'SUN_FLUX' : 'Solar_Radiation',
                             'TERR_MSL' : 'Terrain_Altitude',
                             'THETA' : 'Potential_Temperature'})

        self.folder = folder
        self.filename = filename

        if cfolder is not None:
            self.cfolder = cfolder
            if os.path.exists(os.path.join(self.cfolder,
                                           self.filename + 'CLIPPED')):
                self.cfilename = self.filename + 'CLIPPED'
                self.cfullpath = os.path.join(self.cfolder, self.cfilename)

        self.trajcolor = 'black'
        self.linewidth = 2
        self.name = self.filename + name

    def set_rainstatus(self, rainy_criterion='Rainfall', check_steps=1,
                       rh_threshold=0.8):
        """
        Determines if ``Trajectory`` is 'rainy', sets ``self.rainstatus``
        accordingly.

        The determination is made by examining the user-given criteria
            within the user-given window of timesteps.

        Keyword Arguments
        -----------------
        rainy_criterion : string
            ['Rainfall'|'relative humidity'|'specific humidity']
            'Rainfall' : set ``self.rainstatus`` to ``True`` if rain
            occurs within the indicated ``check_steps``
            'relative humidity' : set ``self.rainstatus`` to ``True`` if
                ``self.relative_humidity`` is above the given ``rh_threshold``
                within ``check_steps``
            'specific humidity' : set ``self.rainstatus`` to ``True`` if
                rain occurs and ``self.specific_humidity`` decreases
        check_steps : integer
            The number of timesteps from the beginning to inspect.
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
            # Test if Rainfall is recorded
            if np.any(self.Rainfall[:check_steps] > 0.0):
                if rainy_criterion == 'Rainfall':
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

    def set_vector(self):
        """
        Calculate mean bearing of ``Trajectory`` and bearings between timesteps

        """

        self.meanvector, self.bearings = ta.tracemean_vector(self.latitude,
                                                             self.longitude)

    def set_relativehumidity(self):
        """
        Calculate relative humidity.  Requires ``Temperature`` and
        ``mixing_ratio`` or ``specific_humidity``

        """

        if not hasattr(self, 'relative_humidity'):
            self.set_mixingratio()
            self.relative_humidity = ta.convert_w2rh(self.mixing_ratio,
                                                     self.Temperature,
                                                     self.pressure)

    def set_mixingratio(self):
        """
        Calculate mixing ratio.  Requires ``specific_humidity`` or
        ``Temperature`` and ``relative_humidity``.

        """

        if not hasattr(self, 'mixing_ratio'):
            if hasattr(self, 'specific_humidity'):
                self.mixing_ratio = ta.convert_q2w(self.specific_humidity)

            elif hasattr(self, 'relative_humidity'):
                self.mixing_ratio = ta.convert_rh2w(self.relative_humidity,
                                                    self.Temperature,
                                                    self.pressure)
            else:
                raise AttributeError('Not enough information to ' +
                                     'calculate mixing ratio!')

    def set_specifichumidity(self):
        """
        Calculate specific humidity.  Requires ``mixing_ratio`` or
        the calculation of ``mixing_ratio`` from ``Temperature`` and
        ``relative_humidity``.

        """

        if not hasattr(self, 'specific_humidity'):
            if hasattr(self, 'mixing_ratio'):
                self.specific_humidity = ta.convert_w2q(self.mixing_ratio)
            elif hasattr(self, 'relative_humidity'):
                self.set_mixingratio()
                self.specific_humidity = ta.convert_w2q(self.mixing_ratio)
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

    def dq_dw_drh(self, qtype):
        """
        Calculate the change in moisture between each timestep

        Works for ``self.specific_humidity``, ``self.mixing_ratio``,
        and ``self.relative_humidity``.

        """
        humidity_dict = {'specific_humidity' : 'dq',
                         'relative_humidity' : 'drh',
                         'mixing_ratio' : 'dw'}

        try:
            humidity = getattr(self, qtype)
        except:
            raise AttributeError('Please calculate ' + qtype +
                                 ' and try again.')

        dhumid_list = []

        for i in range(0, self.sim_length):
            dhumid = humidity[i] - humidity[i + 1]
            dhumid_list.append(dhumid)

        dhumid_arr = np.asarray(dhumid_list).astype(np.float64)
        dhumid_arr = np.pad(dhumid_arr, (0, 1), 'constant',
                            constant_values=(-999.0, -999.0))
        dhumid_arr = np.ma.masked_less_equal(dhumid_arr, -999.0)

        attr_name = humidity_dict[qtype]

        setattr(self, attr_name, dhumid_arr)

    def calculate_moistureflux(self, qtype='specific_humidity'):
        """
        Calculate the moisture flux between each timestep
        using ``self.distance`` and the chosen humidity attribute

        Keyword Arguments
        -----------------
        qtype : string
            ['specific_humidity'|'mixing_ratio']
            Default 'specific_humidity'.

        """

        if not hasattr(self, 'Distance'):
            self.set_distance()

        try:
            humidity = getattr(self, qtype)
        except:
            raise AttributeError('Please calculate ' + qtype +
                                 ' and try again.')

        speed = self.distance / 3600
        mf = speed[1:] * humidity[:-1]

        mf = np.pad(mf, (0, 1), 'constant', constant_values=(-999.0, -999.0))
        self.moistureflux = np.ma.masked_less_equal(mf, -999.0)

    # def last_moistureuptake(self, q_type='specific_humidity'):
    #     """
    #     Method for finding moisture sources according to Gustafsson 2010.

    #     There are three basic scenarios for finding moisure sources.
    #     In the first, the trajectory follows the case outlined in
    #     Gustafsson 2010, where the trajectory experiences a decrease
    #     in q in the target region- the source region is easily
    #     identifiable.  The second is where the parcel shows an increase
    #     in q to the target area
    #     """

    def moistureuptake(self, rainout_threshold, evap_threshold,
                       uptake_window=6, window_overlap=0,
                       vertical_criterion='pbl', pressure_threshold=900.0,
                       mixdepth_factor=1, q_type='specific_humidity'):

        """
        Calculate the moisture flux between each timestep
        using ``self.distance`` and the chosen humidity parameter

        Parameters
        ----------
        rainout_threshold : float
            The level below which changes in humidity are considered to
            indicate loss due to precipitation.  Recommended -0.2 or lower
        evap_threshold : float
            The level above which changes in humidity are considered to
            indicate moisture uptake.  Recommended 0.2 to 0.5
        uptake_window : int
            Default 6.  The length in hours of the period between q0 and qx.
            Assumes that either evaporation or precipitation dominate over
            a short period of time, such as 6 hours.
        window_overlap : int
            Default 0.  The number of hours each ``uptake_window`` should
            overlap the previous one.  A value of 0 indicates that the windows
            should be discrete.  The latest timepoint in the earliest window
            will still start at endpoint + ``uptake_window`` hours.
        vertical_criterion : string
            Default 'pbl'.  ['pbl'|'prs'|'both']
            Criterion for differentiating surficial and other sources of
            moisture.
        pressure_threhold : float
            Default 900.0. The pressure level defined as equivalent
            to the planetary boundary layer (parcel is considered in
            communication with sample)
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

        try:
            initial_timept.extend([self.pressure[-1], self.mixdepth[-1],
                                   self.altitude[-1],
                                   getattr(self, q_type)[-1]])
        except:
            initial_timept.extend([self.pressure[-1], -999.0,
                                   self.altitude[-1],
                                   getattr(self, q_type)[-1]])

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
            try:
                mixdepth = (np.mean(self.mixdepth[ts:ts + uptake_window])
                            * mixdepth_factor)
            except:
                mixdepth = -999.0
                vertical_criterion = 'prs'

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

    def load_reversetraj(self, reverse_dir, fname_end):
        """
        Acquires data from reverse trajectory.

        Parameters
        ----------
        reverse_dir : string
            Location of reverse trajectory

        """

        # Check for directory
        if not os.path.isdir(reverse_dir):
            raise OSError('Reverse trajectory directory does not exist!')

        # Construct filename of reverse trajectory
        self.reversepath = os.path.join(reverse_dir, self.filename + fname_end)

        # Check that this trajectory has reverse trajectory
        if not os.path.exists(self.reversepath):
            raise OSError('File not found: ' + self.reversepath)

        rdatlist, _, _ = hh.load_hysplitfile(self.reversepath)

        # Check that there is only one trajectory in the file
        if len(rdatlist) > 1:
            print 'Multiple-trajectory files not supported!'
            self.reversepath = None
            self.rdata = None
        else:
            self.rdata = rdatlist[0]
            self.rlatitude = self.rdata[:, self.header.index('Latitude')]
            self.rlongitude = self.rdata[:, self.header.index('Longitude')]
            self.raltitude = self.rdata[:,
                                        self.header.index('Altitude')]
            self.rdistance = ta.distance_overearth(self.rlatitude,
                                                   self.rlongitude)
            self.rtotal_dist = ta.sum_distance(self.rdistance)

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
        if not hasattr(self, 'rtotal_dist'):
            raise AttributeError('Reverse trajectory data missing!')
        if not hasattr(self, 'total_distance'):
            self.set_distance()

        # Acquire total horizontal and vertical distances
        r_distance = self.rtotal_dist[-1]
        b_distance = self.total_distance[-1]
        falt_range = np.max(self.raltitude) - np.min(self.raltitude)
        balt_range = np.max(self.altitude) - np.min(self.altitude)

        # Gather coordinates of back trajectory launch location
        # according to both back and forward trajectories.
        site_lats = [self.latitude[0], self.rlatitude[-1]]
        site_lons = [self.longitude[0], self.rlongitude[-1]]

        site_lats = np.asarray(site_lats).astype(np.float64)
        site_lons = np.asarray(site_lons).astype(np.float64)

        # Horizontal distance between back traj launch and forward traj end pts
        site_distance = ta.distance_overearth(site_lats, site_lons)[1]

        # Vertical distance between back traj launch and forward traj end pts
        z_distance = self.altitude[0] - self.raltitude[-1]

        # Distance between points divided by total h or v distance
        self.integ_error_xy = (((site_distance / (r_distance + b_distance))
                               * 100) / 2)
        self.integ_error_z = (((z_distance / (falt_range + balt_range))
                              * 100) / 2)

    def map_traj_scatter(self, cavemap, variable, sizevar=None, **kwargs):
        """
        Scatter plot of ``Trajectory`` data

        Parameters
        ----------
        cavemap : ``Basemap`` instance
            Any ``Basemap`` instance.  For easy map creation, see ``MapDesign``
            class.
        variable : string
            The attribute to plot as a color change
        sizevar : string
            Default ``None``.  The attribute to plot as a marker size change
        **kwargs
            passed to ``traj_scatter()``, then ``cavemap.scatter()`` and
            ``Axes.scatter()``

        Returns
        -------
        cm : ``matplotlib PathCollection`` instance
            Mappable for use in creating colorbars.  Colorbars may be created
            in ``PySPLIT`` using ``make_cbar()`` or ``make_cax_cbar()``

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
        Line plot of ``Trajectory`` path.

        Parameters
        ----------
        cavemap : ``Basemap`` instance
            Any ``Basemap`` instance.  For easy map creation, see ``MapDesign``
            class.
        color : string, tuple
            Default ``None``.  If ``None``, then ``self.color`` will be used.
            Any ``matplotlib``-approved color accepted
        lw : int
            Default ``None``.  If ``None``, then ``self.linewidth``
            will be used.
        **kwargs
            passed to ``traj_plot()`` and then ``cavemap.plot()``,
            and then ``Axes.plot()``

        """

        if color is None:
            color = self.trajcolor
        if lw is None:
            lw = self.linewidth

        mm.traj_path(cavemap, self.longitude, self.latitude,
                     color, lw, **kwargs)
