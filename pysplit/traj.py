from __future__ import division, print_function

# External imports
import os
import numpy as np
import pandas as pd
import geopandas as gp
from shapely.geometry import Point, LineString

# Relative imports within the package
from .hyfile_handler import load_hysplitfile
from .hypath import HyPath


class Trajectory(HyPath):
    """
    Class for processing individual HYSPLIT back trajectories.

    :subclass: of ``HyPath``.

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

        self.data['Temperature_C'] = self.data['Temperature'] - 273.15
        if self.data.get('Mixing_Depth') is None:
            self.data['Mixing_Depth'] = None

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
                sat_vapor = 6.11*(10.0**((7.5*self.data['Temperature_C']) /
                                         (237.7 + self.data['Temperature_C'])))

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
                self.data.calculate_distance()

            self.data['Moisture_Flux'] = None
            self.data.loc[self.index[1]:, 'Moisture_Flux'] = (
                (self.data['Distance_ptp'] / 3600).iloc[1:] *
                self.data.get(humidity).iloc[:-1])


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
        windows = self.data.index[::-interval]

        self.uptake = gp.GeoDataFrame(data=np.empty((windows.size, 13)),
                                      columns=['DateTime', 'Timestep',
                                               'Cumulative_Dist',
                                               'Avg_Pressure',
                                               'Avg_MixDepth', 'q',
                                               'dq_initial', 'dq', 'above',
                                               'below', 'unknown_total',
                                               'above_total', 'below_total'],
                                      dtype=np.float64)

        # Get all the timesteps, set as index
        self.uptake.loc[:, 'Timestep'] = windows[::-1]
        # print('Timestep\n', self.uptake['Timestep'])
        self.uptake.set_index('Timestep', inplace=True, drop=False)

        # fill up with NaNs
        self.uptake.loc[:, ['dq_initial', 'dq', 'above', 'below']] = None

        # (fake) midpoint
        mdpt = interval // 2
        # print('Midpoint', mdpt)

        # Average over the whole window
        for w in windows:
            # print('interval', self.loc[w: w - (interval - 1), 'Timestep'])
            (self.uptake.loc[w, 'Avg_Pressure'],
             self.uptake.loc[w, 'Avg_MixDepth']) = (
                self.data.loc[w: w - (interval - 1),
                              ['Pressure', 'Mixing_Depth']].mean())

        # First timestep
        (self.uptake.loc[windows[0], 'DateTime'],
         self.uptake.loc[windows[0], 'Cumulative_Dist'],
         self.uptake.loc[windows[0], 'q']) = (
            self.data.loc[windows[0], ['DateTime', 'Cumulative_Dist',
                                       humidity]])

        self.uptake.loc[windows[0], 'unknown_total'] = 1.0
        self.uptake.loc[windows[0], ['above_total', 'below_total']] = 0.0

        points.append(self.data.loc[windows[0], 'geometry'])

        for w in windows[1:]:
            (self.uptake.loc[w, 'DateTime'],
             self.uptake.loc[w, 'Cumulative_Dist']) = (
                 self.data.loc[w - mdpt, ['DateTime', 'Cumulative_Dist']])

            self.uptake.loc[w, 'q'] = self.data.loc[w, humidity]

            z = np.mean([pt.z for pt in self.data.loc[w:w - (interval - 1),
                                                      'geometry']])

            points.append(Point([self.data.loc[w - mdpt, 'geometry'].x,
                                 self.data.loc[w - mdpt, 'geometry'].y, z]))

            # set dq initial for timepoints after the earliest:
            self.uptake.loc[w, 'dq_initial'] = (
                self.uptake.loc[w, 'q'] -
                self.uptake.loc[w - interval, 'q'])

        # Set geometry for new gdf
        self.uptake['geometry'] = points

        # Check that mixing depth data actually exists for this trajectory
        if self.data.loc[:, 'Mixing_Depth'].all(None):
            vlim = 'prs'
        else:
            self.uptake.loc[:, 'Avg_MixDepth'] = (
                self.uptake.loc[:, 'Avg_MixDepth'] * mixdepth_factor)

        # Save these results rather than build anew each loop
        is_above = pd.Series([False] * len(self.uptake),
                             index=self.uptake.index, dtype=bool)
        is_below = pd.Series([False] * len(self.uptake),
                             index=self.uptake.index, dtype=bool)

        for wnum, w in enumerate(windows[1:]):
            # wnum is actually the ind of the previous window
            is_above.loc[windows[wnum]] = (
                self.uptake.loc[windows[wnum], 'above'].notnull())
            is_below.loc[windows[wnum]] = (
                self.uptake.loc[windows[wnum], 'below'].notnull())

            if self.uptake.loc[w, 'dq_initial'] > evaporation:
                # set dq
                self.uptake.loc[w, 'dq'] = self.uptake.loc[w, 'dq_initial']

                # Adjust previous fractions
                self.uptake.loc[is_above, 'above'] = (
                    self.uptake.loc[is_above, 'dq'] / self.uptake.loc[w, 'q'])

                self.uptake.loc[is_below, 'below'] = (
                    self.uptake.loc[is_below, 'dq'] / self.uptake.loc[w, 'q'])

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

                fracname_dict = {True: ('below', 'above_total', is_above,
                                        'below_total', is_below),
                                 False: ('above', 'below_total', is_below,
                                         'above_total', is_above)}

                fracs = fracname_dict[is_surface]

                # set new f (is_surface) or e
                self.uptake.loc[w, fracs[0]] = (self.uptake.loc[w, 'dq'] /
                                                self.uptake.loc[w, 'q'])
                # Set new e_total (is_surface) or f_total
                self.uptake.loc[w, fracs[1]] = (
                    self.uptake.loc[fracs[2], 'dq'].sum() /
                    self.uptake.loc[w, 'q'])
                # set new f_total (is_surface) or e_total
                self.uptake.loc[w, fracs[3]] = (
                    (self.uptake.loc[fracs[4], 'dq'].sum() +
                     self.uptake.loc[w, 'dq']) / self.uptake.loc[w, 'q'])
                # Set new unknown fraction
                self.uptake.loc[w, 'unknown_total'] = 1.0 - (
                    self.uptake.loc[w, 'above_total'] +
                    self.uptake.loc[w, 'below_total'])

            else:
                # copy previous total fractions
                self.uptake.loc[w, 'above_total'] = (
                    self.uptake.loc[w - interval, 'above_total'])
                self.uptake.loc[w, 'below_total'] = (
                    self.uptake.loc[w - interval, 'below_total'])
                self.uptake.loc[w, 'unknown_total'] = (
                    self.uptake.loc[w - interval, 'unknown_total'])

                if self.uptake.loc[w, 'dq_initial'] < precipitation:
                    # Adjust previous fractions
                    self.uptake.loc[is_below, 'dq'] = (
                        self.uptake.loc[is_below, 'below'] *
                        self.uptake.loc[w, 'q'])

                    self.uptake.loc[is_above, 'dq'] = (
                        self.uptake.loc[is_above, 'above'] *
                        self.uptake.loc[w, 'q'])

    def load_reversetraj(self, reverse_dir='default', fname_end='REVERSE'):
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

        """
        if not hasattr(self, 'path_r'):

            if reverse_dir is 'default':
                reverse_dir = os.path.join(self.folder, 'reversetraj')

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

            _, path, _, _, multitraj = load_hysplitfile(self.rfullpath)

            if multitraj:
                self.path_r = LineString(
                    [Point(path[self.parcel_ind][i, :]) for i in
                     range(path[self.parcel_ind].shape[0])])
            else:
                self.path_r = LineString(
                    [Point(path[i, :]) for i in range(path.shape[0])])

            # Calculate distance!
            self.calculate_distance(reverse=True)

    def calculate_integrationerr(self):
        """
        Estimate integration error.

        Integration error based on distance between origin and
        reverse trajectory endpoint and the total travel distance.
        Error is in percent.

        """
        if self.data.get('Distance_ptp_r') is None:
            raise AttributeError('Reverse trajectory must be loaded first!')

        if self.data.get('Distance_ptp') is None:
            self.calculate_distance()

        site_distance = self.distance_between2pts(self.path.coords[0],
                                                  self.path_r.coords[-1],
                                                  in_xy=True)

        travel_distance = self.data.loc[:,['Cumulative_Dist',
                                           'Cumulative_Dist_r']].iloc[-1].sum()

        self.integration_error = ((site_distance / travel_distance) * 100) / 2

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
