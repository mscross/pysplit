import numpy as np
import math
import os
import fnmatch
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
import matplotlib.cm as cm
import matplotlib.colors as clrs
from mpl_toolkits.basemap import Basemap
import mapmaker as mm
import hyfile_handler as hh
import traj_accessory as ta

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
        _, self.filename = os.path.split(self.fullpath)

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



    def set_rainstatus(self, rainy_criterion='rainfall', check_steps=1,
                       rh_threshold=0.8):
        """
        Determines if a trajectory is 'rainy', sets self.rainstatus accordingly.

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
                    qx = self.specific_humidity[check_steps*2]
                    dq = q0 - qx

                    if dq < 0:
                        is_rainy = True

        self.rainstatus = is_rainy


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
            self.relative_humidity = self.data[:, self.header.index('RELHUMID')]
        else:
            if hasattr(self, 'mixing_ratio'):
                self.relative_humidity = ta.convert_w2rh(self)
            else:
                self.set_mixingratio()
                self.relative_humidity = ta.convert_w2rh(self)


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
            self.specific_humidity = self.data[:, self.header.index('SPCHUMID')]
        else:
            if hasattr(self, 'mixing_ratio'):
                self.specific_humidity = ta.convert_w2q(self.mixing_ratio)
            elif 'H2OMIXRA' in self.header or hasattr(self,'relative_humidity'):
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

        i = 0

        dq_list = []
        dw_list = []
        drh_list = []

        while i < self.sim_length:
            dq = self.specific_humidity[i]-self.specific_humidity[i+1]
            dw = self.mixing_ratio[i]-self.mixing_ratio[i+1]
            drh = self.relative_humidity[i]=self.relative_humidity[i+1]

            dq_list.append(dq)
            dw_list.append(dw)
            drh_list.append(drh)
            i +=1

        self.dq = np.asarray(dq_list).astype(np.float64)
        self.dw = np.asarray(dw_list).astype(np.float64)
        self.drh = np.asarray(drh_list).astype(np.float64)

        self.dq = np.pad(self.dq, (0,1), 'constant',
                         constant_values=(-999.0, -999.0))
        self.dw = np.pad(self.dw, (0,1), 'constant',
                         constant_values=(-999.0, -999.0))
        self.drh = np.pad(self.drh, (0,1), 'constant',
                         constant_values=(-999.0, -999.0))

        self.dq = np.ma.masked_less_equal(self.dq, -999.0)
        self.dw = np.ma.masked_less_equal(self.dw, -999.0)
        self.drh = np.ma.masked_less_equal(self.drh, -999.0)


    def calculate_moistureflux(self, qtype= 'specific_humidity'):
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

        mf = np.pad(mf, (0,1), 'constant', constant_values=(-999.0, -999.0))
        self.moistureflux = np.ma.masked_less_equal(mf, -999.0)


    def last_moistureuptake(self, q_type='specific_humidity'):
        """
        Method for finding moisture sources according to Gustafsson 2010.

        There are three basic scenarios for finding moisure sources.
            In the first, the trajectory follows the case outlined in
            Gustafsson 2010, where the trajectory experiences a decrease
            in q in the target region- the source region is easily identifiable.
            The second is where the parcel shows an increase in q to
            the target area
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

        for item in moisture_header[:moisture_header.index('total_distance')+1]:
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
        timesteps = range(0, earliest-(uptake_window-1),
                          uptake_window-window_overlap)[::-1]
        print timesteps

        # Cycle through each chunk of time
        for ts in timesteps:

            timepoint = []

            # Get year, month, day, hour, and timestep at the latest
            # point in ts
            for item in moisture_header[:moisture_header.index('timesteps')+1]:
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
                                        moisture_header.index('total_distance')+1]:
                timepoint.append(getattr(self, item)[ts+midpoint])

            # Get the average pressure, mixdepth, and altitude over the window
            pressure = np.mean(self.pressure[ts:ts+uptake_window])
            mixdepth = (np.mean(self.mixdepth[ts:ts+uptake_window])
                        * mixdepth_factor)
            altitude = np.mean(self.altitude[ts:ts+uptake_window])

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
                if (altitude - mixdepth < 0) and (pressure > pressure_threshold):
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
                        updated_frac = mostrecent_dq/q
                        moisture_array[i, e_index] = updated_frac

                    elif moisture_array[i, f_index] > -999.0:
                        mostrecent_dq = moisture_array[i, dq_index]
                        updated_frac = mostrecent_dq/q
                        moisture_array[i, f_index] = updated_frac

                # Initialize current values of e and f
                # Find previous values of dq where there is e and f
                # Get new e_total, f_total
                if below_vert_criteria:
                    e =  -999.0
                    e_inds =  np.nonzero(moisture_array[:, e_index] > -999.0)
                    e_total =  np.sum(moisture_array[e_inds, dq_index])/q

                    f =  dq/q
                    f_inds =  np.nonzero(moisture_array[:, f_index] > -999.0)
                    f_total = (np.sum(moisture_array[f_inds, dq_index])+dq)/q

                else:
                    e =  dq/q
                    e_inds =  np.nonzero(moisture_array[:, e_index]>-999.0)
                    e_total = (np.sum(moisture_array[e_inds, dq_index])+dq)/q

                    f =  -999.0
                    f_inds =  np.nonzero(moisture_array[:, f_index]>-999.0)
                    f_total =  np.sum(moisture_array[f_inds, dq_index])/q

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

            moisture_array = np.concatenate((timepoint, moisture_array), axis=0)

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
        if len(fdatlist) >1:
            print 'Multiple-trajectory files not supported!'
            self.forwardpath = None
            self.fdata = None
        else:
            self.fdata = fdatlist[0]
            self.flatitude = self.fdata[:, self.header.index('Latitude')]
            self.flongitude = self.fdata[:, self.header.index('Longitude')]
            self.fdistance = ta.distance_overearth(self.flatitude,
                                                   self.flongitude)
            self.ftotal_dist = ta.sum_distance(self.fdistance)


    def integration_error(self):
        """
        Estimates the integration error based on distance between back/forward
            start/end points and the total distance traveled.
        Integration error is a minor source of trajectory error.  Resolution
            error should also be checked.

        """

        # Check for necessary data
        if not hasattr(self, 'ftotal_dist'):
            raise AttributeError('Forward trajectory data missing!')
        if not hasattr(self, 'total_distance'):
            self.set_distance()

        f_distance = self.ftotal_dist[-1]
        b_distance = self.total_distance[-1]

        site_lats = [self.latitude[0], self.flatitude[-1]]
        site_lons = [self.longitude[0], self.flongitude[-1]]

        site_lats = np.asarray(site_lats).astype(np.float64)
        site_lons = np.asarray(site_lons).astype(np.float64)

        site_distance = ta.distance_overearth(site_lats, site_lons)[1]

        self.integ_error_est = (site_distance/(f_distance+b_distance))*100



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
        self.directory,_ = os.path.split(traj_object_list[0].fullpath)

        self.map_prefs = {'mapcorners' : [None, None, None, None],
                          'projection' : 'cyl',
                          'standard_pm' : [None, None, None, None],
                          'mapcolor' : 'dark',
                          'shapefiles' : (None, None, None),
                          'labels' : (None, None)}


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


    def get_map_prefs(self):
        """
        Prints the current dictionary of map_prefs

        """

        print 'Current Map Settings:'

        for key, num in zip(self.map_prefs, range(1, 7)):
            print '\t', num, '. ', key, ' : ', self.map_prefs[key]

        print '\n'


    def set_map_prefs(self):
        """
        Adjust the basemap settings.

        Basemap settings
        -----------------
        mapcorners : list of floats
            Used to construct the map view for conic, cylindrical projections.
            Adjusts orientation of map in other projections.
            Lower left longitude, latitude; upper right longitude, latitude.

        projection : string
            Lambert Conformal Conic : 'cconic'
            Albers Equal Area Conic : 'eaconic'
            Equal Area Cylindrical : 'eacyl'
            Equidistant Cylindrical (default) : 'cyl'
            Orthographic : 'ortho'
            Polar Stereographic (conformal) : 'npstere', 'spstere'
            Polar Azimuthal (equal area) : 'nplaea', 'splaea'

        standard_pm : list of floats
            lon_0
                The longitude of central point in all projections except polar
                In polar projections, the longitude that will be oriented
                    to the 12:00 position.
            lat_0
                The latitude of central point in all projections except polar
                In polar projections, this is the boundinglat, which indicates
                    the lowest latitude to appear on the map
            lat_1, lat_2
                Standard parallels in cylindrical and conic projections
                Not necessary to define orthographic or polar projections

        mapcolor : string
            'dark' or 'light'
            The grayscheme of the basemap

        shapefiles : tuple of strings
            Default (None, None, None).
            Shapefile path, color, linewidth.  Leave as deafault to not
                plot shapefiles.

        labels : tuple of strings
            Default (None, None).
            Label group, full or relative path to label file.
            Label group choices:  ['all'|'important'|'justcave'|'cave']
            If path to label file does not exist, the user will be presented
                with several options, including one to make
                and edit a label file

        """

        self.get_map_prefs()

        keep_editing = True
        default_labels = (None, None)
        default_shapefiles = (None, None, None)

        while keep_editing:

            if self.map_prefs['mapcorners'][0] is None:
                to_edit = ['1', '2', '3', '4', '5', '6']

            else:
                print ('Enter the number(s) of the map preference settings ',
                       'you would like to edit: \t')

                to_edit = raw_input()
                to_edit = to_edit.replace(',', ' ').split()

            if '1' in to_edit:
                print 'Enter the following coordinate pairs:\n'
                lowerleft = raw_input('The lat, lon of lower left corner: \t')
                upperright = raw_input('The lat, lon of upper right corner: \t')

                lowerleftlat, lowerleftlon = decompose(lowerleft, 2)
                upperrightlat, upperrightlon = decompose(upperright, 2)

                self.map_prefs['mapcorners'] = [lowerleftlon, lowerleftlat,
                                                upperrightlon, upperrightlat]


            if '2' in to_edit:
                # Edit the projection
                proj_dict = {'1' : 'cconic',
                             '2' : 'eaconic',
                             '3' : 'eacyl',
                             '4' : 'cyl',
                             '5' : 'ortho',
                             '6' : {'N' : 'npstere', 'S' : 'spstere'},
                             '7' : {'N' : 'nplaea',  'S' : 'splaea'}}

                print ('\nChoose a projection from the list by entering ',
                       'its number:\n')
                print ('\t 1. Lambert Conformal Conic\n'
                       '\t 2. Albers Equal Area Conic\n'
                       '\t 3. Equal Area Cylindrical\n'
                       '\t 4. Equidistant Cylindrical (Basemap default)\n'
                       '\t 5. Orthographic\n'
                       '\t 6. Polar Stereographic (Conformal)\n'
                       '\t 7. Polar Azimuthal (Equal-Area)\n')
                proj = raw_input()

                if proj == '6' or proj == '7':
                    print 'You have chosen a polar projection: '
                    hem = raw_input('Please indicate which pole: \t')
                    hem = hem[0].upper()

                    self.map_prefs['projection'] = proj_dict[proj][hem]

                else:
                    self.map_prefs['projection'] = proj_dict[proj]


            if '3' in to_edit:
                # Input the standard parallels and meridians
                print 'Set the orientation of the map: '
                if self.map_prefs['projection'][1] != 'p':

                    lonlat0 = raw_input('Orthographic, Conic, or Cylindrical: '+
                                        ' Provide the central lat, lon: \t')

                    lat_0, lon_0 = decompose(lonlat0, 2)

                    if proj != 'ortho':
                        standard_parallels = raw_input('Conic or Cylindrical: '+
                                                       'Input two standard ' +
                                                       ' parallels: \t')

                        lat_1, lat_2 = decompose(standard_parallels, 2)
                    else:
                        lat_1 = 0.0
                        lat_2 = 0.0
                else:
                    lon_0 = float(raw_input('Polar:  Provide the longitude ' +
                                            'to be oriented vertically at ' +
                                            ' the top of the map: \t'))
                    lat_0 = float(raw_input('\t Provide the bounding latitude: \t'))

                    lat_1 = 0.0
                    lat_2 = 0.0

                self.map_prefs['standard_pm'] = [lon_0, lat_0, lat_1, lat_2]

            if '4' in to_edit:
                mapcolor_opts = {'1' : 'dark',
                                 '2' : 'light'}

                print '\nChoose color scheme for map by entering its number:\n'
                print ('\t 1. dark\n'
                       '\t 2. light\n')

                mc = raw_input()

                self.map_prefs['mapcolor'] = mapcolor_opts[mc]

            if '5' in to_edit:
                yesno_shapefile = raw_input('Plot a shapefile?  Y/N: \t')

                if yesno_shapefile[0].upper() == 'Y':
                    shape_path = raw_input('\nEnter the full or relative ' +
                                            'path to the shapefile: \t')
                    shape_color = raw_input('Enter the shapefile outline ' +
                                            'color (color name or ' +
                                            'hexidecimal color code): \t')
                    shape_lw = int(raw_input('Enter the shapefile ' +
                                             'outline linewidth: \t'))

                    self.map_prefs['shapefiles'] = (shape_path, shape_color,
                                                    shape_lw)
                else:
                    self.map_prefs['shapefiles'] = default_shapefiles

            if '6' in to_edit:
                use_labels = raw_input('Label map?  Y/N: \t')

                if use_labels[0].upper() == 'Y':
                    label_path = raw_input('\nEnter the full or relative ' +
                                           'path to existing or new ' +
                                           'label file: \t')

                    _, label_path = mm.labelfile_reader(label_path)

                    label_dict = {'1' : 'all',
                                  '2' : 'important',
                                  '3' : 'cave',
                                  '4' : 'justcave'}

                    print '\n Available label groups:\n'

                    print ('\t 1. All (countries, oceans, seas, caves)\n'
                           '\t 2. Important (countries, oceans)\n'
                           '\t 3. Cave (countries, oceans, caves)\n'
                           '\t 4. Justcave (caves) \n')

                    lg = raw_input()

                    label_group = label_dict[lg]

                    self.map_prefs['labels'] = (label_group, label_path)

                else:
                    self.map_prefs['labels'] = default_labels


            self.get_map_prefs()

            prefs_ok = raw_input('\nAre map settings ok?  Y/N\t')

            if prefs_ok[0].upper() == 'Y':
                keep_editing = False



    def get_basemap(self, figsize=(15,15), ax=None):
        """
        Takes the current map_prefs and creates the basemap.

        """

        if self.map_prefs['mapcorners'][0] is None:
            self.set_map_prefs()


        cavemap = mm.make_basemap(self.map_prefs['mapcorners'],
                                  self.map_prefs['standard_pm'],
                                  self.map_prefs['projection'], figsize,
                                  self.map_prefs['labels'],
                                  self.map_prefs['shapefiles'],
                                  self.map_prefs['mapcolor'], ax)

        return cavemap


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
            The unique `sort_bytime`s.  Returned if `sort_bytime` is not 'none'.
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


    def grid_trajvar(self, variable, cell_value, grid_res= 0.5,
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
        use_wherebin = Boolean
            Default True.  Use wherebin to assemble a new self.grid, if
            gridding has already occurred

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
        if (not hasattr(self, 'wherebin') or use_wherebin == False
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
            newgrid = np.zeros(self.grid.shape, dtype= self.grid.dtype)

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
            newgrid = np.ma.masked_less_equal(newgrid,-999.0)

            # Set attribute
            self.grid = newgrid

        if normalize:
            gridmax = np.max(self.grid)
            gridmin = np.min(self.grid)

            self.grid = (self.grid-gridmin)/(gridmax-gridmin)


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
                e_rows = np.nonzero(moisture_data[:, e_col]>-999.0)
                mlons = moisture_data[e_rows, lon_col]
                mlats = moisture_data[e_rows, lat_col]
                mdata = moisture_data[e_rows, data_col]

            if 'f' in data_rows:
                f_rows = np.nonzero(moisture_data[:, f_col]>-999.0)
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

            self.mgrid = (self.mgrid-gridmin)/(gridmax-gridmin)


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
                raincount +=1

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
            output = str(traj.fullpath)
            output = output.replace ('\\', '/')
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

        clusteringdata, traj_inds, totalclusters = hh.load_clusterfile(cfile)

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


    def map_data(self, variable, point_size,
                 colorbar_label, add_traj=False, add_scatter=True,
                 traj_color='black', linewidth=2.0, color_max=None,
                 color_min=0.0, colormap='blues', scatter_alpha=1.0,
                 tickdiv=5, figsize=(20,20), variable_transform='linear',
                 sizevar=None, sizevar_transform='linear'):
        """
        Maps the trajectories in trajectory group as scatter/lines.

        Data may be scatter plotted as
        color change and is transformed or not according to user preferences.
        Scatter plot may have second variable plotted as size.  User may choose
        to view trajectories as lines rather than/in addition to scatter.
        Labels, projection, and shapefiles are passed to mm.make_basemap().

        Parameters
        ----------
        variable : string
            Must be in `header`.  Th variable to plot as color change
        point_size : int
            The size of scatter points.  If vary_size == True, then
            point_size is multiplied with the sizevar
        colorbar_label : string
            Label for colorbar

        Keyword Arguments
        -----------------
        add_traj : Boolean
            Default False, trajectory lines not plotted
        add_scatter : Boolean
            Default True, scatter points are plotted
        traj_color : string or tuple
            Default 'black'.
        linewidth : float
            Default 2.0. Width in points of trajectory line.
        color_max : int
            Default None.  The maximum value on the colorbar.
            If None, color_max will be the maximum value of the data.
        color_min : int
            Default 0.0.  The minimum value on the colorbar.
        colormap : string
            Default 'blues'.  ['jet'|'blues'|'anomaly'|'heat'|'earth']
            Passed to a dictionary which retrieves the indicated colormap
        scatter_alpha : float
            Default 1.0.  The opacity of the scatter points
        tickdiv : int
            Default 5. Number of nice ticks on colorbar.  4,5,6 recommended
        figsize : tuple of ints
            Map size
        variable_transform : string
            Default 'linear'.  ['linear'|'sqrt'|'square'|'ln'|'log']
            Determines how data of variable is transformed, if at all
        sizevar : string
            Default None
            Indicates which variable to plot as a change in marker size.
            Only activated if vary_size is True.
        sizevar_transform : string
            Default 'linear'.  ['linear'|'sqrt'|'square'|'ln'|'log']
            Determines how data of sizevar is transformed, if at all

        Returns
        -------
        cavemap : plot
            plot of HYSPLIT trajectory data

        """

        # Make basemap from self.map_prefs
        cavemap = self.get_basemap(figsize=figsize)

        if add_traj:

            colorlist = []

            if type(traj_color) is tuple:
                min_height = traj_color[0]
                max_height = traj_color[1]

            for traj in self.trajectories:
                # Calculate the correct grayscale value of the trajectory
                # based on given minimum and maximum height
                # Lower altitude trajectories are darker gray

                if type(traj_color) is tuple:
                    grey_val = ((((traj.altitude_t0 - min_height)/
                                  (max_height   - min_height))*(1.0-0.25))+0.25)

                    colorlist.append(str(grey_val))

                # If (min, max) not given, trajectories will be all traj_color
                else:
                    colorlist.append(traj_color)

            for traj, color in zip(self.trajectories, colorlist):
                cavemap.plot(traj.longitude, traj.latitude, color=color,
                             linewidth=linewidth, latlon=True, zorder=19)

        if add_scatter:

            data, lons, lats = data_prep(self, variable, variable_transform)
            colormap = mm.get_colormap(colormap)

            print data

            if color_max == None:
                color_max = data.max()

            if sizevar is not None:
                data2 = data_prep(self, sizevar, sizevar_transform)
                data2 = data2 * point_size
                point_size = data2

            cm = cavemap.scatter(lons, lats, c=data, s=point_size,
                                 cmap=colormap, vmin=color_min, vmax=color_max,
                                 latlon=True, zorder=20, edgecolor='none',
                                 alpha = scatter_alpha)


            cbar = cavemap.colorbar(cm, location='bottom', pad='5%')

            cbar.locator = tk.MaxNLocator(tickdiv, integer=False)
            cbar.ax.tick_params(labelsize=16)
            cbar.update_ticks()

            cbar.set_label(colorbar_label, labelpad=10, fontsize=16)

        return cavemap


    def map_data_verticalplot(self, variable, point_size, colorbar_label,
                              verticallevels, vertplot_xlabel, vertplot_ylabel,
                              vertplot_yticks, hspace, add_scatter='all',
                              add_traj='none', traj_color='black',
                              color_max=None, color_min=0.0,
                              vertplot_timelim=(-168,0), vertplot_ylim=(None),
                              colormap='blues',
                              cax_pos=[0.925, 0.15, 0.03, 0.65],
                              scatter_alpha=1.0, tickdiv=5, figsize=(20,20),
                              relative_mapsize=8, variable_transform='linear',
                              vary_size=False, sizevar=None,
                              sizevar_transform='linear'):

        """
        Similar to map_data(), except below the map, a user-specified variable
            is plotted in time v vertical coordinate space.  The colorbar is
            placed to the left of the two plots.

        Adjustment of colorbar axis location/size and hspace likely necessary
            to achieve the best layout of map, vertical plot, and colorbar.

        Parameters
        ----------
        variable : string
            The variable to plot as color change
        point_size : int
            The size to make each scatter point.  If vary_size, then
            point_size is multiplied with the sizevar
        colorbar_label : string
            Label for colorbar
        verticallevels : string
            ['altitude'|'pressure']  The vertical coordinates
        vertplot_xlabel : string
            X-axis label for the vertical v time plot
        vertplot_ylabel : string
            Y-axis label  for the vertical v time plot
        vertplot_yticks : tuple of ints
            The (major, minor) tick multiples on vertical axis of the vertplot
        hspace : float
            Spacing adjustment between map and vertplot.  Recommended
            to start at -0.45 and tweak from there

        Keyword Arguments
        -----------------
        add_traj : string
            ['none'|'map'|'plot'|'all']
            Default 'none'.  Add trajectory lines to map and/or vertplot,
            or neither.
        add_scatter : string
            ['none'|'map'|'plot'|'all']
            Default 'all'.  Add scatter points to map and/or vertplot,
            or neither
        traj_color : string or tuple
            If string, all trajectories will be plotted the same color
            If tuple (minimum altitude, maximum altitude), then trajectories
            will be plotted in grayscale by altitude (lower = darker)
        color_max : int
            Default None.  The emaximum value on the colorbar.
            If None, color_max will be the maximum value of the data.
        color_min : int
            Default 0.0.  The minimum value on the colorbar.
        vertplot_timelim
        vertplot_ylim
            Default 0.0.  The minimum value on the colorbar.
        colormap : string
            Default 'blues'.  ['jet'|'blues'|'anomaly'|'heat'|'earth']
            Passed to a dictionary which retrieves the indicated colormap
        cax_pos : list of floats
            Default [0.925, 0.15, 0.03, 0.65].  The position of the colorbar.
            [left, bottom, width, height]
        scatter_alpha : float
            Default 1.0 (opaque).  Alpha value of the scatter points.
        tickdiv : int
            Default 5. Number of nice ticks on the colorbar.  4,5,6 recommended
        figsize : tuple of ints
            Size of figure
        relative_mapsize : int

        variable_transform : string
            Default 'linear'.  ['linear'|'sqrt'|'square'|'ln'|'log']
            Determines how data of variable is transformed, if at all
        vary_size : Boolean
            Default False.  (True|False)
            If True, indicates to vary size as well as color of scatter points
        sizevar : string
            Default None
            Indicates which variable to plot as a change in marker size.
            Only activated if vary_size is True.
        sizevar_transform : string
            Default 'linear'.  ['linear'|'sqrt'|'square'|'ln'|'log']
            Determines how data of sizevar is transformed, if at all

        """

        # Initialize figure and relative sizes of plots
        fig = plt.figure(figsize=figsize)
        rows = relative_mapsize
        mapspan = rows - 1

        ax0 = plt.subplot2grid((rows, 1), (0,0), rowspan=mapspan)
        ax1 = plt.subplot2grid((rows, 1), (mapspan, 0))

        # Invert axis if plotting on pressure levels
        if verticallevels == 'pressure':
            print 'YUP'
            ax1.invert_yaxis()

        # Initialize the basemap on the first axis
        cavemap = self.get_basemap(figsize=figsize, ax=ax0)

        # Prepare to draw individual trajectories
        if add_traj != 'none':
            # Initialize lists of x,y, and color arrays
            colorlist = []
            vertlist = []

            for traj in self.trajectories:

                # Calculate the correct grayscale value of the trajectory
                # based on given minimum and maximum height
                # Lower altitude trajectories are darker gray
                if type(traj_color) is tuple:
                    min_height = traj_color[0]
                    max_height = traj_color[1]

                    grey_val = ((((traj.altitude_t0 - min_height)/
                                  (max_height - min_height))*(1.0-0.25)))

                    colorlist.append(str(grey_val))

                # If (min, max) not given, trajectories will be all the given
                else:
                    colorlist.append(traj_color)

            for traj, color in zip(self.trajectories, colorlist):
                # Plot trajectories on map
                if add_traj != 'plot':
                    cavemap.plot(traj.longitude, traj.latitude, color=color,
                                 latlon=True, zorder=19)

                # Plot trajectories on vertplot
                if add_traj != 'map':
                    ax1.plot(traj.timesteps, getattr(traj, verticallevels),
                             color=color, zorder=19)

        vert = []
        traveltime = []

        for traj in self.trajectories:
            v = getattr(traj, verticallevels)
            t = getattr(traj, 'timesteps')

            v = v.tolist()
            t = t.tolist()

            vert.extend(v)
            traveltime.extend(t)


        if add_scatter != 'none':

            # Prepare primary (color) data, get the colormap
            data, lons, lats = data_prep(self, variable, variable_transform)
            colormap = mm.get_colormap(colormap)

            if color_max is None:
                color_max = data.max()

            if vary_size:
                data2 = data_prep(self, sizevar, sizevar_transform)
                data2 = data2 * point_size
                point_size = data2

            if add_scatter != 'plot':
                # Plot on map
                cm = cavemap.scatter(lons, lats, c=data, s=point_size,
                                     cmap=colormap, vmin=color_min,
                                     vmax=color_max, latlon=True, zorder=20,
                                     edgecolor='none', alpha=scatter_alpha)

            if add_scatter != 'map':
                cm = ax1.scatter(traveltime, vert, c=data, s=point_size,
                                 cmap=colormap, edgecolor='none',
                                 vmin=color_min, vmax=color_max,
                                 alpha=scatter_alpha, zorder=16)

            # Initialize colorbar
            cax = fig.add_axes([cax_pos[0], cax_pos[1],
                                cax_pos[2], cax_pos[3]])
            cbar = plt.colorbar(cm, cax=cax)

            # Format colorbar ticks and labels
            cbar.locator = tk.MaxNLocator(tickdiv, integer=False)
            cbar.ax.tick_params(labelsize=16)
            cbar.update_ticks()
            cbar.set_label(colorbar_label, labelpad=10, fontsize=16,
                           rotation=270, verticalalignment='bottom')

        # Set x and y limits on vertplot
        ax1.set_xlim(vertplot_timelim[1], vertplot_timelim[0])

        if vertplot_ylim==None:
            ax1.set_ylim(0, np.max(vert))
            if verticallevels == 'pressure':
                ax1.invert_yaxis()
        else:
            ax1.set_ylim(vertplot_ylim[0], vertplot_ylim[1])

        # Format y axis ticks and all ticklabel sizes for vertplot
        ymajorlocator = tk.MultipleLocator(vertplot_yticks[0])
        yminorlocator = tk.MultipleLocator(vertplot_yticks[1])

        ax1.yaxis.set_major_locator(ymajorlocator)
        ax1.yaxis.set_minor_locator(yminorlocator)
        ax1.tick_params(labelsize=16)

        ax1.set_xlabel(vertplot_xlabel, fontsize=16)
        ax1.set_ylabel(vertplot_ylabel, fontsize=16)

        plt.subplots_adjust(hspace=hspace)

        return cavemap


    def map_moisture(self, uptake, scale, point_size,
                     colorbar_label, color_max=None, color_min=None,
                     scatter_alpha=1.0, tickdiv=5, figsize=(20,20)):
        """
        Plot moisture uptakes as a scatter plot.

        Plot all points considered for moisture uptake; plot all uptakes in
            fractional or absolute values, plot uptakes below or uptakes
            above vertical criteria in fractional or absolute values

        Parameters
        ----------
        uptake : string
            ['both'|'above'|'below'|'all points']
        scale : string
            ['absolute dq'|'absolute dqi|'fractional']
        point_size : int
            The size to make each scatter point.
        colorbar_label : string
            Label for colorbar

        Keyword Arguments
        -----------------
        color_max : string
            Default None.  Maximum value on colorbar.  If not given and
            `scale` is absolute, value is introspected from data.  If not given
            and `scale` is fractional, color_max is 1.0.
        color_max : string
            Default None.  Minimum value on colorbar.  If not given and
            `scale` is absolute, value is introspected from data.  If not given
            and `scale` is fractional, color_min is 0.0.
        scatter_alpha : float
            Default 1.0 (opaque).  Alpha value of the scatter points.
        tickdiv : int
            Default 5. Number of nice ticks on the colorbar.  4,5,6 recommended
        figsize : tuple of ints
            Size of figure

        Returns
        -------
        cavemap

        """

        # Create Basemap from self.map_prefs
        cavemap = self.get_basemap(figsize=figsize)

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

        if color_max is None:
            if 'absolute' in scale:
                color_max=0.0
                for tr in self.trajectories:
                    datamax = np.max(tr.masked_sources[:, data_col])
                    if datamax > color_max:
                        color_max = datamax
            else:
                color_max = 1.0

        if color_min is None:
            if 'absolute' in scale:
                color_min=0.0
                for tr in self.trajectories:
                    datamin = np.min(tr.masked_sources[:, data_col])
                    if datamin < color_min:
                        color_min = datamin
            else:
                color_min = 0.0

        for tr in self.trajectories:
            if data_rows is None:
                lons = tr.masked_sources[:, lon_col]
                lats = tr.masked_sources[:, lat_col]
                data = tr.masked_sources[:, data_col]
                cm = cavemap.scatter(lons, lats, c=data, s=point_size,
                                     cmap=cmap,vmin=color_min, vmax=color_max,
                                     latlon=True, zorder=20, edgecolor='none')
            else:
                if 'e' in data_rows:
                    e_rows = np.nonzero(tr.masked_sources[:, e_col]>-999.0)
                    lons = tr.masked_sources[e_rows, lon_col]
                    lats = tr.masked_sources[e_rows, lat_col]
                    data = tr.masked_sources[e_rows, data_col]
                    cm = cavemap.scatter(lons, lats, c=data, s=point_size,
                                         cmap=cmap,vmin=color_min,
                                         vmax=color_max, latlon=True,
                                         zorder=20, edgecolor='none')
                if 'f' in data_rows:
                    f_rows = np.nonzero(tr.masked_sources[:, f_col]>-999.0)
                    lons = tr.masked_sources[f_rows, lon_col]
                    lats = tr.masked_sources[f_rows, lat_col]
                    data = tr.masked_sources[f_rows, data_col]
                    cm = cavemap.scatter(lons, lats, c=data, s=point_size,
                                         cmap=cmap,vmin=color_min,
                                         vmax=color_max, latlon=True,
                                         zorder=20, edgecolor='none')

        cbar = cavemap.colorbar(cm, location='bottom', pad='5%')

        cbar.locator = tk.MaxNLocator(tickdiv, integer=False)
        cbar.ax.tick_params(labelsize=16)
        cbar.update_ticks()

        cbar.set_label(colorbar_label, labelpad=10, fontsize=16)

        return cavemap


    def gridmap(self, colorbar_label, usecontourf=False, mapcount=False,
                ismoisture=False, color_max=None, color_min=None,
                colormap='blues', scatter_alpha=1.0, tickdiv=5,
                figsize=(20,20)):
        """
        Create a pcolor plot of gridded data.

        Parameters
        ----------
        colorbar_label : string
            Label for colorbar

        Keyword Arguments
        -----------------
        ismoisture : Boolean
            Access gridded moisture uptake data if True.  If False,
            access regular gridded data.
        color_max : string
            Default None.  Maximum value on colorbar.  If not given and
            `scale` is absolute, value is introspected from data.  If not given
            and `scale` is fractional, color_max is 1.0.
        color_min : string
            Default None.  Minimum value on colorbar.  If not given and
            `scale` is absolute, value is introspected from data.  If not given
            and `scale` is fractional, color_min is 0.0.
        colormap : string
            Default 'blues'.  ['jet'|'blues'|'anomaly'|'heat'|'earth']
            Passed to a dictionary which retrieves the indicated colormap
        scatter_alpha : float
            Default 1.0 (opaque).  Alpha value of the scatter points.
        tickdiv : int
            Default 5. Number of nice ticks on the colorbar.  4,5,6 recommended
        figsize : tuple of ints
            Size of figure

        Returns
        -------
        cavemap

        """

        # Make Basemap from self.map_prefs
        cavemap = self.get_basemap(figsize=figsize)

        colormap = mm.get_colormap(colormap)

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

        if usecontourf:
            cm = cavemap.contourf(x,y, data,
                                  np.linspace(color_min, color_max, num=11),
                                  cmap=colormap, latlon=True, zorder=20)
        else:
            cm = cavemap.pcolormesh(x, y, data, latlon=True, cmap=colormap,
                                    vmin=color_min, vmax=color_max, zorder=20)

        cbar = cavemap.colorbar(cm, location='bottom', pad='5%')

        cbar.locator = tk.MaxNLocator(tickdiv, integer=False)
        cbar.ax.tick_params(labelsize=16)
        cbar.update_ticks()

        cbar.set_label(colorbar_label, labelpad=10, fontsize=16)

        return cavemap



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
        # filename = 'C'+str(self.clusternumber) +'_'+str(total_clusters)+'mean.tdump'
        # endpoints_file = os.path.join(endpoints_dir, filename)

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



class ClusterGroup:
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


    def map_clusters(self, color_var='mean_mf', color_min=0.0,
                     color_max=None, color_transform='linear',
                     color_adjust=1.0, width_var='count',
                     width_transform='linear', width_adjust=1.0,
                     colormap='blues', figsize=(20,20)):
        """
        Plots mean cluster paths with color and width scaled to some variable.

        Creates a map where each cluster trajectory is plotted with a color
            determined by the mean variable value and a width by the trajectory
            count.  Color and width scaled according to user preferences.

        Keyword Arguments
        -----------------
        color_var : string
            Default 'mean_mf'.  ['total_raint0'|'total_mft1'|'mean_q'|
            'mean_w'|'mean_mf'|'mean_mft1'|'mean_raint0'|'total_mf']
        color_min : float
            Default 0.0.  The minimum value of the color variable.
            Used in normalization of color variable data to 0, 1 scale
        color_max : float
            Default None.  The maximum value of the color variable.
            If None, color_max is the data max.  Used in normalization of color
            variable data to 0, 1 scale.
        color_transform : string
            Default 'linear'.  ['linear'|'sqrt'|'square'|'ln'|'log']
            Determines how color data is transformed, if at all
        width_var : string
            Default 'count'.  ['absolute_count'|'relative_count'|
            'total_raint0'|'total_mft1'|'mean_q'|
            'mean_w'|'mean_mf'|'mean_mft1'|'mean_raint0'|'total_mf']
        width_transform : string
            Default 'linear'.  ['linear'|'sqrt'|'square'|'ln'|'log']
            Determines how width data is transformed, if at all
        width_adjust : float
            Default 1.0.  Used to adjust linewidths to a reasonable value
            while maintaining relationships between linesizes.
        colormap : string
            Default 'blues'.  ['jet'|'blues'|'anomaly'|'heat'|'earth']
            Passed to a dictionary which retrieves the indicated colormap
        figsize : tuple of ints
            Size of figure

        Returns
        -------
        cavemap : plot
            Map of cluster paths drawn to indicate number of component
            trajectories and average of some data

        """

        # Make basemap from self.map_prefs
        cavemap = self.get_basemap(figsize=figsize)

        # Acquire transforms
        c_transf = get_transform(color_transform)
        w_transf = get_transform(width_transform)

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
                    tmp.append((w/self.totaltrajcount)*100)
                widthvar_list = tmp
        else:
            for cluster in self.clusters:
                if not hasattr(cluster, width_var):
                    cluster.set_meanvar()
                widthvar_list.append(getattr(cluster, width_var))

        # Transform and adjust widths, transform colors
        i = 0
        while i < self.totalclusters:
            widthvar_list[i] = (w_transf[0](widthvar_list[i]))*width_adjust
            colorvar_list[i] = c_transf[0](colorvar_list[i])
            i = i + 1

        # Get color max
        if color_max is None:
            color_max = max(colorvar_list)

        # Set up mapping of scalar data to RGBA values from given colormap
        colormap = mm.get_colormap(colormap)
        cnorm = clrs.Normalize(vmin=color_min, vmax=color_max)
        scalarmap = cm.ScalarMappable(norm=cnorm, cmap=colormap)

        # Obtain index array of linewidths so thick lines plot last
        plot_order = np.argsort(widthvar_list, kind='mergesort')

        for i in plot_order:
            # Map values to RGBA values from given colormap, use resulting color
            color = scalarmap.to_rgba(colorvar_list[i])
            cavemap.plot(self.clusters[i].longitude, self.clusters[i].latitude,
                         linewidth=widthvar_list[i], color=color, zorder=19,
                         latlon=True)

        return cavemap


    def trajplot_clusterplot(self, colors=None, orientation='horizontal',
                             cluster_lw='relative', traj_lw=3, figsize=(20,20)):

        """
        Plots trajectories on one map and cluster paths on another.

        Trajectories belonging to a cluster will be plotted as lines of
            the same color, and in the cluster plot the cluster mean path
            will also be the same color.  The linewidth of the cluster mean
            path will reflect the number of member trajectories unless
            specified otherwise.

        Keyword Arguments
        -----------------
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
        traj_lw : int
            Default 3.  The linewdith of each trajectory.

        Returns
        -------
        fig : figure instance
            Figure with two subplots
        trajmap : Basemap instance
            Left/top subplot: map with trajectory paths plotted on it.
            Trajectory colors correspond to cluster
        clusmap : Basemap instance
            Right/bottom subplot: map with cluster mean paths plotted on it.
            Cluster color corresponds to that of member trajectories.
            Linewidth is given value or the number of member trajectories.

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
            color_tmp=[]
            for c in colors:
                color_tmp.append(c[0])
            colors = color_tmp

        # Set cluster mean path width by member trajectory count if not given
        if cluster_lw is 'relative':
            print 'relative'
            cluster_lw = []
            for cluster in self.clusters:
                cluster_lw.append((cluster.trajcount/float(self.totaltrajcount))*100.0)

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
        trajmap = self.get_basemap(figsize=figsize, ax=ax_t)
        clusmap = self.get_basemap(figsize=figsize, ax=ax_c)

        # Plot
        for cluster, color, lw in zip(self.clusters, colors, cluster_lw):
            clusmap.plot(cluster.longitude, cluster.latitude, color=color,
                         linewidth=lw, latlon=True, zorder=19, ax=ax_c)

            for traj in cluster.trajectories:
                trajmap.plot(traj.longitude, traj.latitude, color=color,
                             linewidth=traj_lw, latlon=True, zorder=19, ax=ax_t)

        return fig, trajmap, clusmap


def data_prep(trajgroup, variable, transform):
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

    # Get the functions necessary for transforming the data
    transforms = get_transform(transform)

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

    data = np.ma.masked_less_equal(datarray,-999.0)
    lons = lonarray
    lats = latarray

    data = transforms[1](data)

    return data, lons, lats


def get_transform(transform):
    """
    Dictionary of transforms.  Keys are strings representing transforms,
        values are functions

    Parameters
    ----------
    transform : string
        How the data is to be transformed.
        ['sqrt'|'linear'|'log'|'ln'|'square']

    Returns
    -------
    transform_dict : list of functions
        The functions necessary to adjust data to specified transform

    """

    # 0 is used for cluster mapping, 1 for regular trajectory mapping
    transform_dict = {'sqrt'  : [math.sqrt, np.sqrt],
                      'linear': [do_nothing, do_nothing],
                      'log'   : [math.log10, np.log10],
                      'ln'    : [math.log, np.log],
                      'square': [square_it, square_it]}

    transforms = transform_dict[transform]

    return transforms


def square_it(x):
    """
    """
    x = x**2
    return x


def do_nothing(x):
    """
    Does nothing
    """
    return x


def decompose(stringinput, desired_num):
    """
    Takes a string of numbers, removes commas and separates into a list,
        then converts items into floats, which may or may not
        be unpacked when decompose() is called.

    Parameters
    ----------
    stringinput : string
    desired_num : int

    Returns
    -------
    float_input : list of floats

    """

    # Split
    stringinput = stringinput.replace(',', ' ').split()

    # Check that we have the number of things that we want
    if len(stringinput) != desired_num:
        raise ValueError('You have entered the wrong number of inputs.' +
                         '  Please enter ' + str(desired_num) + ' items.')

    # Convert strings to floats
    float_input = []
    for item in stringinput:
        float_input.append(float(item))

    return float_input