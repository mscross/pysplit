from __future__ import division, print_function

# External imports
import os
import numpy as np
import pandas as pd
import geopandas as gp
from shapely.geometry import Point, LineString
from subprocess import call
import fnmatch

# Relative imports within the package
from .hyfile_handler import load_hysplitfile
from .hypath import HyPath
from .trajectory_generator import (_populate_control, _try_to_remove,
                                   _day2filenum, _cliptraj)


class Trajectory(HyPath):
    """
    Class for processing individual HYSPLIT back trajectories.

    :subclass: of ``HyPath``.

    """

    def __init__(self, trajdata, pathdata, datetime, trajheader, folder,
                 filename, cfolder, multitraj):
        """
        Initialize ``Trajectory``.

        Parameters
        ----------
        trajdata : (M, N) ndarray of floats
            The data array corresponding to the along-path data of a
            single HYSPLIT trajectory.
        pathdata : (M, 3) ndarray of floats
            The longitude, latitude, and altitude information corresponding
            to ``trajdata``.
        datetime : (M) pandas.DatetimeIndex
            The dates and times corresponding to each timestep of
            ``trajdata`` and ``pathdata``.
        trajheader : list of N strings
            The column headers for ``trajdata``.
        folder : string
            Path information of HYSPLIT file.
        filename : string
            Path information of HYSPLIT file.
        cfolder : string
            Location of corresponding clipped HYSPLIT file, if it exists.
        multitraj : Boolean
            If True, is from a file containing multiple trajectories
        """
        HyPath.__init__(self, trajdata, pathdata, datetime, trajheader)

        self.data.rename(columns={'AIR_TEMP': 'Temperature',
                                  'PRESSURE': 'Pressure',
                                  'RAINFALL': 'Rainfall',
                                  'MIXDEPTH': 'Mixing_Depth',
                                  'RELHUMID': 'Relative_Humidity',
                                  'H2OMIXRA': 'Mixing_Ratio',
                                  'SPCHUMID': 'Specific_Humidity',
                                  'SUN_FLUX': 'Solar_Radiation',
                                  'TERR_MSL': 'Terrain_Altitude',
                                  'THETA': 'Potential_Temperature'},
                         inplace=True)

        # Not everyone has Temperature output
        try:
            self.data['Temperature_C'] = self.data['Temperature'] - 273.15
        except KeyError:
            self.data['Temperature_C'] = None
            self.data['Temperature'] = None

        self.folder = folder
        self.filename = filename
        self.fullpath = os.path.join(self.folder, self.filename)

        self.trajid = self.fullpath + str(trajdata[0, 0])

        self.trajcolor = 'black'
        self.linewidth = 2
        self.multitraj = multitraj
        self.parcel_num = trajdata[0, 0]
        # False until proven otherwise

    def __hash__(self):
        """Magic hash method."""
        return hash(self.trajid)

    def __eq__(self, other):
        """Magic eq method."""
        if isinstance(other, self.__class__):
            return self.trajid == other.trajid
        return NotImplemented

    def set_rainstatus(self, rainy_criterion='Rainfall', check_steps=0,
                       threshold=0.0):
        """
        Determine if ``Trajectory`` produced rain during indicated timesteps.

        Doesn't look at change in specific humidity.

        Parameters
        ----------
        rainy_criterion : string
            Default 'Rainfall'.  The attribute to check to determine if rainy.
        check_steps : int or slice
            Default 0.  The timestep(s) (the index, not positional arguments)
            to check for rainy-ness.
        threshold : float
            Default 0.0.  Minimum value for rainyness.  0.8 suggested for
            ``Relative_Humidity``.

        """
        if self.data.get(rainy_criterion) is None:
            raise KeyError(rainy_criterion, " not in Trajectory.data.columns")

        self.rainy = False
        if np.any(self.data.loc[check_steps, rainy_criterion] > threshold):
            self.rainy = True

    def calculate_rh(self):
        """
        Calculate ``Relative_Humidity`` from ``Mixing_Ratio``.

        Will not execute if ``Relative_Humidity`` is already present.
        Calculation ultimately requires either ``Mixing_Ratio`` or
        ``Specific_Humidity`` as original along-path output variables.

        """
        # Check for existence of relative humidity and mixing ratio
        if self.data.get('Relative_Humidity') is None:
            if self.data.get('Mixing_Ratio') is None:
                raise KeyError('Calculate mixing ratio first!')
            else:
                # Convert mixing ratio to relative humidity
                sat_vapor = 6.11 * (10.0**((7.5 * self.data['Temperature_C']) /
                                           (237.7 +
                                            self.data['Temperature_C'])))

                sat_w = 621.97 * (sat_vapor / (self.data['Pressure'] -
                                               sat_vapor))

                self.data['Relative_Humidity'] = ((self.data['Mixing_Ratio'] /
                                                   sat_w) * 100.0)

    def calculate_w(self, calc_using):
        """
        Calculate ``Mixing_Ratio``.

        Use either ``Relative_Humidity`` or ``Specific_Humidity``.
        Will not execute if ``Mixing_Ratio`` is already present.

        Parameters
        ----------
        calc_using : string
            The humidity parameter to use to calculate ``Mixing_Ratio``.
            [``Relative_Humidity``|``Specific_Humidity``]

        """
        if self.data.get('Mixing_Ratio') is None:
            if self.data.get(calc_using) is None:
                raise KeyError(calc_using, ' does not exist.')
            else:
                func_dict = {'Relative_Humidity': self._convert_rh2w,
                             'Specific_Humidity': self._convert_q2w}

                func_dict[calc_using]()

    def calculate_sh(self):
        """
        Calculate ``Specific_Humidity`` from ``Mixing_Ratio``.

        Will not execute if ``Specific_Humidity`` is already present.
        Calculation ultimately requires either ``Mixing_Ratio`` or
        ``Relative_Humidity`` as original along-path output variables.

        """
        if self.data.get('Specific_Humidity') is None:
            if self.data.get('Mixing_Ratio') is None:
                raise KeyError('Calculate mixing ratio first!')
            else:
                w_kg = self.data['Mixing_Ratio'] / 1000
                self.data['Specific_Humidity'] = (w_kg / (w_kg + 1)) * 1000

    def calculate_moistureflux(self, humidity='Specific_Humidity'):
        """
        Calculate ``Moisture_Flux``.

        Moisture flux between each timestep, uses ``Distance_ptp``
        and the indicated humidity type (``humidity``).

        Parameters
        ----------
        humidity :  string
            Default 'Specific_Humidity'. The humidity parameter used to
            calculate ``Moisture_Flux``.
            [``Relative_Humidity``|``Specific_Humidity``]

        """
        if self.data.get(humidity) is None:
            print('Calculate ', humidity, ' first!')
        else:
            if self.data.get('Distance_ptp') is None:
                self.calculate_distance()

            self.data['Moisture_Flux'] = None
            self.data.loc[self.data.index[1]:, 'Moisture_Flux'] = (
                (self.data['Distance_ptp'] / 3600).iloc[1:] *
                self.data.get(humidity).iloc[:-1])

    def moisture_uptake(self, precipitation, evaporation, interval=6,
                        starting_timepoint=None, vlim='prs',
                        pressure_level=900.0, mixdepth_factor=1,
                        humidity='Specific_Humidity'):
        """
        Moisture uptakes for back trajectories.

        Parameters
        ----------
        precipitation : float
            Suggested -0.2.  The change in humidity below (not inclusive) which
            precipitation is considered to have occurred.
        evaporation : float
            Suggested 0.2 to 0.5.  The change in humidity above
            (not inclusive) which evaporation is considered to have occurred
        interval : int
            Default 6. The length of a calculation window in timesteps.
            -HYSPLIT by default has timesteps 1 hr apart.
            -Assumes that either evaporation or precipitation dominate over
            a short period of time (like 6 hours).
            -Windows are counted UP from the earliest timepoint (or given
            `starting_timepoint`).  If the length in timesteps of the
            trajectory is not equally divisible by `interval`, then to end at
            0 you must change `starting_timepoint`.

            Example:  a 30-hour back trajectory with a default
                `starting_timepoint` will set its initial conditions to the
                conditions at -30.  Then uptakes/decreases will be calculated
                for the -29 to -24, -23 to -18, -17 to -12, -11 to -6,
                and -5 to 0 windows by finding the average pressure, altitude
                for each window, and comparing the humidity at -24 to -30,
                -18 to -24, -12 to -18, -6 to -12, and 0 to -6.
                The `uptake` `GeoDataFrame` will have
                information at -30, -24, -18, -12, -6, and 0.
        starting_timepoint : int
            Default None.  The timepoint at which to start the calculation.
            If None, then calculation will start at the end of the trajectory.
        vlim : string
            Default 'pbl'.  ['pbl'|'prs'|'both']
            Criterion for distinguishing surficial and other moisture sources:
            below the planetary boundary layer, below a given
            ``pressure_level``, or both.  Pressure and altitude are averaged
            over the entire ``interval`` (e.g., from -29 to -24)
        pressure_level : float
            Default 900.0.  The pressure level defined as equivalent to the
            planetary boundary layer.
        mixdepth_factor : int or float
            Default 1.  The value by which to adjust the parcel average
            mixed layer depth.  Use if mixed layer depth seems to be under-
            or over-estimated.
        humidity : string
            Default 'Specific_Humidity'.  The humidity parameter to use to
            calculate moisture uptakes.

        """
        # Check for variables
        if self.data.get(humidity) is None:
            raise ValueError(self.trajid, ' does not have ', humidity)

        use_mixdepth = True
        if self.data.get('Mixing_Depth') is None:
            use_mixdepth = False
            vlim = 'prs'

        # Init and set up uptake GeoDataFrame with initial and meteo/geom
        windows = self._get_uptake_windows(starting_timepoint, interval)

        self._init_uptake_gdf(windows, use_mixdepth)
        self._setup_uptake_md_prs(windows, interval, use_mixdepth,
                                  mixdepth_factor)

        mdpt = interval // 2

        self._setup_uptake(windows, interval, humidity, mdpt)
        self._setup_uptake_geom(windows, interval, mdpt)

        boollist = [False] * len(self.uptake)
        is_above = pd.Series(boollist, index=self.uptake.index, dtype=bool)
        is_below = pd.Series(boollist, index=self.uptake.index, dtype=bool)

        checkerdict = {'prs': self._surfcheck_prs, 'pbl': self._surfcheck_pbl,
                       'both': self._surfcheck_both}

        surface_checker = checkerdict[vlim]

        for wnum, w in enumerate(windows[1:]):
            # wnum is the ind of the previous window
            i = windows[:wnum + 1]
            is_above.loc[i] = self.uptake.loc[i, 'above'].notnull()
            is_below.loc[i] = self.uptake.loc[i, 'below'].notnull()

            if self.uptake.loc[w, 'dq_initial'] > evaporation:
                # set dq
                self.uptake.loc[w, 'dq'] = self.uptake.loc[w, 'dq_initial']
                self._adjust_previous_fracs(is_above, is_below, w)
                is_surface = surface_checker(w, pressure_level)
                self._set_new_fracs(is_above, is_below, w, is_surface)

            else:
                # No fractional change, copy previous total fractions
                for t in ['above_total', 'below_total', 'unknown_total']:
                    self.uptake.loc[w, t] = self.uptake.loc[w - interval, t]

                if self.uptake.loc[w, 'dq_initial'] < precipitation:
                    self._adjust_previous_dq(is_above, is_below, w)

    def load_clippedtraj_data(self, clipped_dir='default',
                              fname_end='CLIPPED'):
        """
        Load folder, filename of clipped version of trajectory.

        Parameters
        ----------
        clipped_dir : string
            The location of the clipped trajectories.  Default is a subfolder
            in ``self.folder`` named 'clippedtraj'
        fname_end : string
            Default 'CLIPPED'. Clipped trajectory filename is ``self.filename``
            + fname_end

        """
        if clipped_dir is 'default':
            clipped_dir = os.path.join(self.folder, 'clippedtraj')

        if not os.path.isdir(clipped_dir):
            raise OSError('Clipped trajectory directory does not exist!')

        self.cfolder = clipped_dir
        cfullpath = os.path.join(self.cfolder,
                                 self.filename + fname_end)
        if os.path.exists(cfullpath):
            self.cfilename = self.filename + fname_end
            self.cfullpath = cfullpath
        else:
            raise OSError(cfullpath, ' not found.')

    def generate_clippedtraj(self, clipped_dir='default'):
        """
        Generate the clipped version of the original trajectory file.

        The clipped trajectory file contains only path and time information.
        Does not require a HYSPLIT installation to create, only the original
        trajectory file.

        Parameters
        ----------
        clipped_dir : string
            Full or relative path to the reverse trajectory directory.
            'default' refers to ``self.folder`` + 'clippedtraj'

        """
        orig_dir = os.getcwd()

        if clipped_dir is 'default':
            clipped_dir = os.path.join(self.folder, 'clippedtraj')

        if not os.path.isdir(clipped_dir):
            os.mkdir(os.path.join(clipped_dir))

        try:
            os.chdir(self.folder)

            _cliptraj(clipped_dir, self.filename)

        finally:
            os.chdir(orig_dir)

    def load_reversetraj(self, reverse_dir='default', fname_end='REVERSE',
                         reload_rtraj=False):
        """
        Load reverse trajectory.

        Load as a LineString, then put distance info into ``self.data``.
        Multi-trajectory files supported (maybe).

        Parameters
        ----------
        reverse_dir : string
            The location of the reverse trajectories.  Default is a subfolder
            in ``self.folder`` named 'reversetraj'.
        fname_end : string
            Default 'REVERSE'. Reverse trajectory filename is ``self.filename``
            + fname_end.  This keyword included to grandfather in trajectories
            calculated with previous ``PySPLIT`` versions
            (``fname_end`` = 'FORWARD').
        reload_rtraj : Boolean
            Default False.  If True, will reload the reverse trajectory.

        """
        if reload_rtraj or not hasattr(self, 'path_r'):

            if reverse_dir is 'default':
                reverse_dir = os.path.join(self.folder, 'reversetraj')
            if not os.path.isdir(reverse_dir):
                raise OSError('Reverse trajectory directory does not exist!')

            self.rfolder = reverse_dir
            rfullpath = os.path.join(self.rfolder, self.filename + fname_end)

            if os.path.exists(rfullpath):
                self.rfilename = self.filename + fname_end
                self.rfullpath = rfullpath
            else:
                raise OSError(rfullpath, ' not found.')

            _, path, _, _, multitraj = load_hysplitfile(self.rfullpath)

            badtraj = False

            if multitraj:
                badlens = []
                badinds = []
                path_r = LineString(
                    [Point(path[self.parcel_num - 1][i, :]) for i in
                     range(path[self.parcel_num - 1].shape[0])])
                for pnum, pr in enumerate(path_r):
                    if len(pr.xy[0]) != len(self.data.index):
                        badtraj = True
                        badlens.append(len(pr.xy[0]) - 1)
                        badinds.append(pnum)

                if badtraj:
                    verb = 'are'
                    if len(badlens) == 1:
                        verb = 'is'
                    args = (self.trajid, badinds, verb, badlens,
                            len(self.data.index) - 1)
                    print('''Trajectory {} has bad reverse trajectories: \n\t
                          {} {} {} hours instead of {} hours'''.format(*args))

            else:
                path_r = LineString(
                    [Point(path[i, :]) for i in range(path.shape[0])])

                if len(path_r.xy[0]) != len(self.data.index):
                    badtraj = True

                    args = (self.trajid, len(path_r.xy[0]) - 1,
                            len(self.data.index) - 1)
                    print('''Trajectory {} has a bad reverse trajectory: \n\t
                          {} hours instead of {} hours'''.format(*args))

            # Calculate distance!
            if badtraj:
                self.rtraj_ok = False
            else:
                self.rtraj_ok = True
                self.path_r = path_r
                self.calculate_distance(reverse=True)

    def generate_reversetraj(self, hysplit_working, meteo_dir,
                             reverse_dir='default',
                             meteo_interval='weekly',
                             hysplit="C:\\hysplit4\\exec\\hyts_std"):
        r"""
        Generate the reverse trajectory.  Requires HYSPLIT installation.

        The reverse trajectory begins at the endpoint of the original
        trajectory, and runs the opposite direction in time.  The two paths
        are ideally indistinguishable when plotted.

        Parameters
        ----------
        hysplit_working : string
            Full or relative path to the HYSPLIT working directory.
        meteo_dir : string
            Full or relative path to the location of the meteorology files.
        reverse_dir : string
            Full or relative path to the reverse trajectory directory.
            Default refers to ``self.folder`` + 'reversetraj'
        meteo_interval : string
            Default 'weekly'.  ['monthly'|semimonthly'|'daily'|'weekly'].
            Whether the meteorlogy files used to calculate trajectory are
            monthly, semi-monthly, weekly, or daily files.
        hysplit : string
            Default "C:\\hysplit4\\exec\\hyts_std".  The location of the
            "hyts_std" executable that generates trajectories.  This is the
            default location for a typical PC installation of HYSPLIT.

        """
        orig_dir = os.getcwd()

        if reverse_dir is 'default':
            reverse_dir = os.path.join(self.folder, 'reversetraj')

        if not os.path.isdir(reverse_dir):
            os.mkdir(os.path.join(reverse_dir))

        reversetrajname = self.filename + 'REVERSE'
        final_rtrajpath = os.path.join(reverse_dir, reversetrajname)

        y = self.data.DateTime.dt.year.iloc[-1]
        m = self.data.DateTime.dt.month.iloc[-1]
        d = self.data.DateTime.dt.day.iloc[-1]
        h = self.data.DateTime.dt.hour.iloc[-1]

        coordinates = (self.data.geometry.iloc[-1].y,
                       self.data.geometry.iloc[-1].x)
        alt = self.data.geometry.iloc[-1].z

        if alt >= 10000.0:
            alt = 9999
        run = self.data.index[-1] * -1

        year_str = '{:02}'.format(int(str(y)[-2:]))

        # Introspect meteorology files
        if not hasattr(self, 'meteorology_files'):
            self.get_meteo_files(meteo_dir, meteo_interval)

        try:
            os.chdir(hysplit_working)

            _try_to_remove('CONTROL')
            _try_to_remove(reversetrajname)
            _try_to_remove(final_rtrajpath)

            _populate_control(coordinates, year_str, m, d, h, alt, meteo_dir,
                              self.meteorology_files, run, 'CONTROL',
                              reversetrajname)

            call(hysplit)

            os.rename(reversetrajname, final_rtrajpath)

        finally:
            os.chdir(orig_dir)

    def get_meteo_files(self, meteo_dir, meteo_interval):
        """
        Get a list of the meteorology files used to calculate ``Trajectory``.

        Results stored as ``self.meteorology_files``

        Parameters
        ----------
        meteo_dir : string
            Full or relative Path to the location of the meteorology files.
        meteo_interval : string
            The time coverage of the meteorological files.
            ['semimonthly'|'daily'|'weekly']

        """
        meteopatterns = []
        meteofiles = []

        orig_dir = os.getcwd()
        meteo_interval = meteo_interval[0].lower()

        mon_dict = {'1': 'jan', '2': 'feb', '3': 'mar', '4': 'apr',
                    '5': 'may', '6': 'jun', '7': 'jul', '8': 'aug',
                    '9': 'sep', '10': 'oct', '11': 'nov', '12': 'dec'}

        with open(self.fullpath, 'r') as trajfile:
            contents = trajfile.readlines()

            for line in contents[1:]:
                if 'OMEGA' in line:
                    break

                parts = line.split()[1:4]

                year = "{:02}".format(int(parts[0]))
                month = mon_dict[parts[1]]
                day = _day2filenum(meteo_interval, parts[2])

                filestring = '*' + month + '*' + year + '*' + day

                meteopatterns.append(filestring)

        try:
            os.chdir(meteo_dir)

            _, _, files = next(os.walk('.'))

            for pattern in meteopatterns:
                for each_file in files:
                    if fnmatch.fnmatch(each_file, pattern):
                        meteofiles.append(each_file)
                        break

        finally:
            os.chdir(orig_dir)

        if len(meteofiles) == 0:
            raise OSError('No meteorology files found.')
        print(len(meteofiles))

        self.meteorology_files = meteofiles

    def calculate_integrationerr(self):
        """
        Estimate integration error.

        Absolute integration error (``self.integration_error_abs``)
        is half the loop closure distance, which is the distance between
        ``Trajectory`` (``self``) origin and reverse trajectory endpoint.
        Relative integration error is the loop closure distance divided by the
        total distance traveled by reverse and original trajectories, divided
        by 2 and reported as a percentage.

        """
        if not hasattr(self, 'rtraj_ok'):
            raise AttributeError('Reverse trajectory must be loaded first!')

        if self.rtraj_ok:

            if self.data.get('Distance_ptp') is None:
                self.calculate_distance()

            site_distance = self.distance_between2pts(self.path.coords[0],
                                                      self.path_r.coords[-1],
                                                      in_xy=True)

            travel_distance = self.data.loc[:,
                ['Cumulative_Dist', 'Cumulative_Dist_r']].iloc[-1].sum()

            self.integration_error = ((site_distance / travel_distance) *
                                      100) / 2
            self.integration_error_abs = site_distance / 2

        else:
            print('''Integration error calculation skipped for
                  Trajectory {}'''.format(self.trajid))

    def _convert_rh2w(self):
        """
        Convert ``Relative_Humidity`` to ``Mixing_Ratio``.

        Only called by ``self.calculate_w()``, private.

        """
        sat_vapor = 6.11 * (10.0 ** ((7.5 * self.data['Temperature_C']) /
                                     (237.7 + self.data['Temperature_C'])))

        sat_w = 621.97 * (sat_vapor / (self.data['Pressure'] - sat_vapor))

        self.data['Mixing_Ratio'] = (
            self.data['Relative_Humidity'] / 100.0) * sat_w

    def _convert_q2w(self):
        """
        Convert ``Specific_Humidity`` to ``Mixing_Ratio``.

        Only called by ``self.calculate_w()``, private.

        """
        q_kg = self.data['Specific_Humidity'] / 1000

        self.data['Mixing_Ratio'] = (q_kg / (1 - q_kg)) * 1000

    def _get_uptake_windows(self, starting_timepoint, interval):
        """Define uptake windows for self.moisture_uptake()."""
        if starting_timepoint is None:
            # Gives -164, -158, -152, -146 ... 0 etc
            windows = self.data.index[::-interval]
        else:
            index = self.data.index.tolist()
            if starting_timepoint > 0:
                raise ValueError('`starting_timepoint` must be negative')
            try:
                i = index.index(starting_timepoint)
            except ValueError:
                print('`starting_timepoint` of ' +
                      str(starting_timepoint) + ' not found in index ' +
                      'of trajectory ' +
                      self.trajid + '. \nDefaulting to end of trajectory at ' +
                      str(self.data.index[-1]))
                windows = self.data.index[::-interval]
            else:
                windows = self.data.index[:i + 1][::-interval]

        return windows

    def _init_uptake_gdf(self, windows, use_mixdepth):
        """Initialize self.uptake with or without mixing depth."""
        if use_mixdepth:
            columns = ['DateTime', 'Timestep', 'Cumulative_Dist',
                       'Avg_Pressure', 'Avg_MixDepth', 'q', 'dq_initial', 'dq',
                       'above', 'below', 'unknown_total', 'above_total',
                       'below_total']
        else:
            columns = ['DateTime', 'Timestep', 'Cumulative_Dist',
                       'Avg_Pressure', 'q', 'dq_initial', 'dq',
                       'above', 'below', 'unknown_total', 'above_total',
                       'below_total']

        self.uptake = gp.GeoDataFrame(data=np.empty((windows.size, 13)),
                                      columns=columns, dtype=np.float64)

        # Get all the timesteps, set as index
        self.uptake.loc[:, 'Timestep'] = windows[::-1]
        self.uptake.set_index('Timestep', inplace=True, drop=False)

        # fill up with NaNs
        self.uptake.loc[:, ['dq_initial', 'dq', 'above', 'below']] = None

    def _setup_uptake_md_prs(self, windows, interval, use_mixdepth,
                             mixdepth_factor):
        """
        Initialize pressure and maybe mixing depth columns in self.uptake.

        debug statement:
        print('interval', self.loc[w: w - (interval - 1), 'Timestep'])

        """
        if use_mixdepth:
            for w in windows:
                (self.uptake.loc[w, 'Avg_Pressure'],
                 self.uptake.loc[w, 'Avg_MixDepth']) = (
                    self.data.loc[w: w - (interval - 1),
                                  ['Pressure', 'Mixing_Depth']].mean())

            self.uptake.loc[:, 'Avg_MixDepth'] = (
                self.uptake.loc[:, 'Avg_MixDepth'] * mixdepth_factor)
        else:
            for w in windows:
                self.uptake.loc[w, 'Avg_Pressure'] = (
                    self.data.loc[w: w - (interval - 1), ['Pressure']].mean())

    def _setup_uptake(self, windows, interval, humidity, mdpt):
        """Initialize datetime, distance, humidity, dqi, q in self.uptake."""
        (self.uptake.loc[windows[0], 'DateTime'],
         self.uptake.loc[windows[0], 'Cumulative_Dist'],
         self.uptake.loc[windows[0], 'q']) = (
            self.data.loc[windows[0], ['DateTime', 'Cumulative_Dist',
                                       humidity]])

        self.uptake.loc[windows[0], 'unknown_total'] = 1.0
        self.uptake.loc[windows[0], ['above_total', 'below_total']] = 0.0

        for w in windows[1:]:
            (self.uptake.loc[w, 'DateTime'],
             self.uptake.loc[w, 'Cumulative_Dist']) = (
                 self.data.loc[w - mdpt, ['DateTime', 'Cumulative_Dist']])

            self.uptake.loc[w, 'q'] = self.data.loc[w, humidity]

            self.uptake.loc[w, 'dq_initial'] = (
                self.uptake.loc[w, 'q'] -
                self.uptake.loc[w - interval, 'q'])

    def _setup_uptake_geom(self, windows, interval, mdpt):
        """Initialize self.uptake.geometry."""
        points = []

        points.append(self.data.loc[windows[0], 'geometry'])

        for w in windows[1:]:
            z = np.mean([pt.z for pt in self.data.loc[w:w - (interval - 1),
                                                      'geometry']])
            points.append(Point([self.data.loc[w - mdpt, 'geometry'].x,
                                 self.data.loc[w - mdpt, 'geometry'].y, z]))

        self.uptake['geometry'] = points[::-1]

    def _surfcheck_prs(self, w, pressure_level):
        """Check if surficial uptake, pressure criterion."""
        if self.uptake.loc[w, 'Avg_Pressure'] > pressure_level:
            is_surface = True
        else:
            is_surface = False

        return is_surface

    def _surfcheck_pbl(self, w, pressure_level):
        """Check if surficial uptake, planetary boundary layer criterion."""
        if (self.uptake.loc[w, 'geometry'].z <
                self.uptake.loc[w, 'Avg_MixDepth']):
            is_surface = True
        else:
            return is_surface

    def _surfcheck_both(self, w, pressure_level):
        """Check if surficial uptake, both vertical criteria."""
        prs = self._surface_checker_prs(w, pressure_level)
        pbl = self._surface_checker_pbl(w, pressure_level)

        if prs and pbl:
            is_surface = True
        else:
            is_surface = False

        return is_surface

    def _adjust_previous_dq(self, is_above, is_below, w):
        """Adjust dq to account for water loss due to precip."""
        self.uptake.loc[is_below, 'dq'] = (
            self.uptake.loc[is_below, 'below'] *
            self.uptake.loc[w, 'q'])

        self.uptake.loc[is_above, 'dq'] = (
            self.uptake.loc[is_above, 'above'] *
            self.uptake.loc[w, 'q'])

    def _adjust_previous_fracs(self, is_above, is_below, w):
        """Adjust previous fracs to account for water gain due to evap."""
        self.uptake.loc[is_above, 'above'] = (
            self.uptake.loc[is_above, 'dq'] / self.uptake.loc[w, 'q'])

        self.uptake.loc[is_below, 'below'] = (
            self.uptake.loc[is_below, 'dq'] / self.uptake.loc[w, 'q'])

    def _set_new_fracs(self, is_above, is_below, w, is_surface):
        """Set current fractions to account for water gain from evap."""
        if is_surface:
            self.uptake.loc[w, 'below'] = (
                self.uptake.loc[w, 'dq'] / self.uptake.loc[w, 'q'])

            self.uptake.loc[w, 'above_total'] = (
                self.uptake.loc[is_above, 'dq'].sum() /
                self.uptake.loc[w, 'q'])

            self.uptake.loc[w, 'below_total'] = (
                (self.uptake.loc[is_below, 'dq'].sum() +
                 self.uptake.loc[w, 'dq']) / self.uptake.loc[w, 'q'])

        else:
            self.uptake.loc[w, 'above'] = (
                self.uptake.loc[w, 'dq'] / self.uptake.loc[w, 'q'])

            self.uptake.loc[w, 'below_total'] = (
                self.uptake.loc[is_below, 'dq'].sum() /
                self.uptake.loc[w, 'q'])

            self.uptake.loc[w, 'above_total'] = (
                (self.uptake.loc[is_above, 'dq'].sum() +
                 self.uptake.loc[w, 'dq']) / self.uptake.loc[w, 'q'])

        self.uptake.loc[w, 'unknown_total'] = 1.0 - (
            self.uptake.loc[w, 'above_total'] +
            self.uptake.loc[w, 'below_total'])
