from __future__ import division, print_function
import os
import numpy as np
import geopandas as gp
from shapely.geometry import Point, LineString
import hyfile_handler as hh
from hypath import HyPath


class Trajectory(HyPath):
    """
    Class for processing individual HYSPLIT back trajectories.
    Subclass of pandas DataFrame.

    """

    def __init__(self, trajdata, pathdata, datetime, trajheader, folder,
                 filename, cfolder):
        """
        Initialize ``Trajectory``.

        Parameters
        ----------
        trajdata : (M, N) ndarray of floats
            The data array corresponding to the along-path data of a
            single HYSPLIT trajectory.
        pathdata : (M, 3) ndarray of floats
            The latitude, longitude, and altitude information corresponding
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

        """

        HyPath.__init__(self, trajdata, pathdata, datetime, trajheader)

        self.rename(columns={'AIR_TEMP': 'Temperature',
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

        self['Temperature_C'] = self['Temperature'] - 273.15
        if self.get('Mixing_Depth') is None:
            self['Mixing_Depth'] = None

        self.folder = folder
        self.filename = filename
        self.fullpath = os.path.join(self.folder, self.filename)

        self.trajid = self.fullpath + str(trajdata[0, 0])

        if cfolder is not None:
            self.cfolder = cfolder
            if os.path.exists(os.path.join(self.cfolder,
                                           self.filename + 'CLIPPED')):
                self.cfilename = self.filename + 'CLIPPED'
                self.cfullpath = os.path.join(self.cfolder, self.cfilename)

        self.trajcolor = 'black'
        self.linewidth = 2

    def __hash__(self):
        return hash(self.trajid)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.trajid == other.trajid
        return NotImplemented

    def set_rainstatus(self, rainy_criterion='Rainfall', check_steps=0,
                       threshold=0.0):
        """
        Determines if ``Trajectory`` produced rain during indicated timesteps.
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
        if self.get(rainy_criterion) is None:
            raise KeyError(rainy_criterion, " is not in Trajectory.columns")

        self.rainy = False

        if np.any(self.loc[check_steps, rainy_criterion] > threshold):
            self.rainy = True

    def calculate_rh(self):
        """
        Calculate ``Relative_Humidity`` from ``Mixing_Ratio``.
        Will not execute if ``Relative_Humidity`` is already present.

        Calculation ultimately requires either ``Mixing_Ratio`` or
        ``Specific_Humidity`` as original along-path output variables.

        """
        # Check for existence of relative humidity and mixing ratio
        if self.get('Relative_Humidity') is None:
            if self.get('Mixing_Ratio') is None:
                raise KeyError('Calculate mixing ratio first!')
            else:
                # Convert mixing ratio to relative humidity
                sat_vapor = 6.11 * (10.0 ** ((7.5 * self['Temperature_C']) /
                                    (237.7 + self['Temperature_C'])))

                sat_w = 621.97 * (sat_vapor / (self['Pressure'] - sat_vapor))

                self['Relative_Humidity'] = ((self['Mixing_Ratio'] /
                                              sat_w) * 100.0)

    def calculate_w(self, calc_using):
        """
        Calculate ``Mixing_Ratio`` using ``Relative_Humidity`` or
        ``Specific_Humidity``.
        Will not execute if ``Mixing_Ratio`` is already present.

        Parameters
        ----------
        calc_using : string
            The humidity parameter to use to calculate ``Mixing_Ratio``.
            [``Relative_Humidity``|``Specific_Humidity``]

        """
        if self.get('Mixing_Ratio') is None:
            if self.get(calc_using) is None:
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
        if self.get('Specific_Humidity') is None:
            if self.get('Mixing_Ratio') is None:
                raise KeyError('Calculate mixing ratio first!')
            else:
                w_kg = self['Mixing_Ratio'] / 1000
                self['Specific_Humidity'] = (w_kg / (w_kg + 1)) * 1000

    def calculate_moistureflux(self, humidity='Specific_Humidity'):
        """
        Calculate ``Moisture_Flux`` between each timestep using
        ``Distance`` and the indicated humidity type (``humidity``).

        Parameters
        ----------
        humidity :  string
            The humidity parameter to use to calculate ``Moisture_Flux``.
            [``Relative_Humidity``|``Specific_Humidity``]

        """

        if self.get(humidity) is None:
            print('Calculate ', humidity, ' first!')
        else:
            if self.get('Distance') is None:
                self.calculate_distance()

            self['Moisture_Flux'] = np.empty((self['Distance'].size))
            self['Moisture_Flux'].iloc[0] = None
            self['Moisture_Flux'].iloc[1:] = ((self['Distance'] /
                                               3600).iloc[1:] *
                                              self.get(humidity).iloc[:-1])

    def moisture_uptake(self, precipitation, evaporation, interval=6,
                        vlim='pbl', pressure_level=900.0,
                        mixdepth_factor=1, humidity='Specific_Humidity'):
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
            Default 6. The number of hourly timesteps between humidity checks.
            Assumes that either evaporation or precipitation dominate over
            a short period of time (like 6 hours).
            Example:  a 30-hour back trajectory will have an initial
                humidity data point at -30, then check again at -24, -18, -12,
                -6, and 0.
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
        points = []

        # Gives 164, 158, 152, 146 ... 0 etc
        windows = self.index[::-interval]

        self.uptake = gp.GeoDataFrame(data=np.empty((windows.size, 13)),
                                      columns=['DateTime', 'Timestep',
                                               'Cumulative_Dist',
                                               'Avg_Pressure',
                                               'Avg_MixDepth', 'q',
                                               'dq_initial', 'dq', 'e', 'f',
                                               'unknown_frac', 'e_total',
                                               'f_total'], dtype=np.float64)

        # Get all the timesteps, set as index
        self.uptake.loc[:, 'Timestep'] = windows[::-1]
        self.uptake.loc.set_index('Timestep', inplace=True)

        # fill up with NaNs
        self.uptake.loc[:, ['dq_initial', 'dq', 'e', 'f']] = None

        # (fake) midpoint
        mdpt = len(windows) // 2

        # Average over the whole window
        for w in windows:
            self.uptake.loc[w, ['Avg_Pressure', 'Avg_MixDepth']] = (
                self.loc[w: w - interval, ['Pressure', 'Mixing_Depth']].mean())

        # First timestep
        self.uptake.loc[windows[0], ['DateTime', 'Cumulative_Dist', 'q']] = (
            self.loc[windows[0], ['DateTime', 'Cumulative_Dist', humidity]])

        self.uptake.loc[windows[0], 'unknown_frac'] = 1.0
        self.uptake.loc[windows[0], ['e_total', 'f_total']] = 0.0

        points.append(self.loc[windows[0], 'geometry'])

        for w in windows[1:]:
            self.uptake.loc[w, ['DateTime', 'Cumulative_Dist', 'q']] = (
                self.loc[w - mdpt, ['DateTime', 'Cumulative_Dist', humidity]])

            z = np.mean([pt.z for pt in self.loc[w:w - interval, 'geometry']])

            points.append(Point([self.loc[w - mdpt, 'geometry'].x,
                                 self.loc[w - mdpt, 'geometry'].y, z]))

            # set dq initial for timepoints after the earliest:
            self.uptake.loc[w, 'dq_initial'] = (
                self.uptake.loc[w, 'q'] -
                self.uptake.loc[w - interval, 'q'])

        # Set geometry for new gdf
        self.uptake['geometry'] = points

        # Check that mixing depth data actually exists for this trajectory
        if self.uptake.loc[:, 'Mixing_Depth'].all(None):
            vlim = 'prs'
        else:
            self.uptake.loc[:, 'Avg_MixDepth'] = (
                self.uptake.loc[:, 'Avg_MixDepth'] * mixdepth_factor)

        for wnum, w in enumerate(windows[1:]):
            is_e = self.uptake.loc[windows[:wnum + 1], 'e'].notnull()
            is_f = self.uptake.loc[windows[:wnum + 1], 'f'].notnull()

            if self.uptake.loc[w, 'dq_initial'] > evaporation:
                # set dq
                self.uptake.loc[w, 'dq'] = self.uptake.loc[w, 'dq_initial']

                # Adjust previous fractions
                self.uptake[is_e].loc[:, 'e'] = (
                    self.uptake[is_e].loc[:, 'dq'] / self.uptake.loc[w, 'q'])

                self.uptake[is_f].loc[:, 'f'] = (
                    self.uptake[is_f].loc[:, 'dq'] / self.uptake.loc[w, 'q'])

                is_surface = False
                if vlim is 'prs':
                    if self.uptake.loc[w, 'Avg_Pressure'] > pressure_level:
                        is_surface = True

                elif vlim is 'pbl':
                    if (self.uptake.loc[w, 'geometry'].z <
                            self.uptake.loc[w, 'Avg_MixDepth']):
                        is_surface = True

                else:
                    if (self.uptake.loc[w, 'Avg_Pressure'] > pressure_level and
                            (self.uptake.loc[w, 'geometry'].z <
                             self.uptake.loc[w, 'Avg_MixDepth'])):
                        is_surface = True

                fracname_dict = {True: ('f', 'e_total', is_e,
                                        'f_total', is_f),
                                 False: ('e', 'f_total', is_f,
                                         'e_total', is_e)}

                fracs = fracname_dict[is_surface]

                # set new f (is_surface) or e
                self.uptake.loc[w, fracs[0]] = (self.uptake.loc[w, 'dq'] /
                                                self.uptake.loc[w, 'q'])
                # Set new e_total (is_surface) or f_total
                self.uptake.loc[w, fracs[1]] = (
                    self.uptake[fracs[2]].loc[:, 'dq'].sum() /
                    self.uptake.loc[w, 'q'])
                # set new f_total (is_surface) or e_total
                self.uptake.loc[w, fracs[3]] = (
                    (self.uptake[fracs[4]].loc[:, 'dq'].sum() +
                     self.uptake.loc[w, 'dq']) / self.uptake.loc[w, 'q'])
                # Set new unknown fraction
                self.uptake.loc[w, 'unknown_frac'] = 1.0 - (
                    self.uptake.loc[w, 'e_total'] +
                    self.uptake.loc[w, 'f_total'])

            else:
                # copy previous total fractions
                self.uptake.loc[w, 'e_total'] = (
                    self.uptake.loc[w - interval, 'e_total'])
                self.uptake.loc[w, 'f_total'] = (
                    self.uptake.loc[w - interval, 'f_total'])
                self.uptake.loc[w, 'unknown_frac'] = (
                    self.uptake.loc[w - interval, 'unknown_frac'])

                if self.uptake.loc[w, 'dq_initial'] < precipitation:
                    # Adjust previous fractions
                    self.uptake[is_f].loc[:, 'dq'] = (
                        self.uptake[is_f].loc[:, 'f'] *
                        self.uptake.loc[w, 'q'])

                    self.uptake[is_e].loc[:, 'dq'] = (
                        self.uptake[is_e].loc[:, 'e'] *
                        self.uptake.loc[w, 'q'])

    def load_reversetraj(self, reverse_dir, fname_end='REVERSE'):
        """
        Loads reverse trajectory as a LineString, then puts distance info
        into geodataframe.

        Multi-trajectory files supported

        """
        if not hasattr(self, 'path_r'):

            if not os.path.isdir(reverse_dir):
                raise OSError('Reverse trajectory directory does not exist!')

            self.rfolder = reverse_dir
            if os.path.exists(os.path.join(self.rfolder,
                                           self.filename + fname_end)):
                self.rfilename = self.filename + fname_end
                self.rfullpath = os.path.join(self.rfolder, self.rfilename)
            else:
                raise OSError('File ', self.filename,
                              fname_end, ' not found.')

            _, path, _, _, multitraj = hh.load_hysplitfile(self.rfullpath)

            if multitraj:
                self.path_r = LineString(
                    [Point(path[self.parcel_ind][i, :]) for i in
                     range(path[self.parcel_ind].shape[0])])
            else:
                self.path_r = LineString(
                    [Point(path[i, :]) for i in range(path.shape[0])])

            # Calculate distance!
            self.calculate_distance(reverse=True)

    def integration_error(self):
        """
        Estimate integration error.

        """
        if self.get('Distance_r') is None:
            raise AttributeError('Reverse trajectory must be loaded first!')

        if self.get('Distance') is None:
            self.calculate_distance()

        site_distance = self.distance_between2pts(self.path.coords[0],
                                                  self.path.coords[-1])

        travel_distance = self.loc[:, ['Cumulative_Dist',
                                       'Cumulative_Dist_r']].iloc[-1].sum()

        self.integration_error = ((site_distance / travel_distance) * 100) / 2

    def _convert_rh2w(self):
        """
        Private function for converting ``Relative_Humidity`` to
        ``Mixing_Ratio``.

        Only called by ``self.calculate_w()``.

        """
        sat_vapor = 6.11 * (10.0 ** ((7.5 * self['Temperature_C']) /
                                     (237.7 + self['Temperature_C'])))

        sat_w = 621.97 * (sat_vapor / (self['Pressure'] - sat_vapor))

        self['Mixing_Ratio'] = (self['Relative_Humidity'] / 100.0) * sat_w

    def _convert_q2w(self):
        """
        Private function for converting ``Specific_Humidity`` to
        ``Mixing_Ratio``.

        Only called by ``self.calculate_w()``.

        """
        q_kg = self['Specific_Humidity'] / 1000

        self['Mixing_Ratio'] = (q_kg / (1 - q_kg)) * 1000
