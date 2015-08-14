import math
import numpy as np
import geopandas as gp
from shapely.geometry import Point, LineString


class HyPath(gp.GeoDataFrame):
    """
    Class for initializing HySPLIT trajectories and cluster paths.

    Subclass of GeoPandas GeoDataFrame.
    Superclass of Trajectory and Cluster

    """

    def __init__(self, alongpath, pathdata, datetime, header):
        """
        Initialize GeoDataFrame and path.

        Parameters
        ----------
        alongpath : (M, N) ndarray of floats
            The data array corresponding to the along-path data of a
            single HYSPLIT trajectory.
        pathdata : (M, 3) ndarray of floats
            The latitude, longitude, and altitude information corresponding
            to ``alongpath``.
        datetime : (M) pandas.DatetimeIndex
            The dates and times corresponding to each timestep of
            ``alongpath`` and ``pathdata``.
        header : list of N strings
            The column headers for ``alongpath``.
        """

        pts = [Point(pathdata[i, :]) for i in range(pathdata.shape[0])]

        gp.GeoDataFrame.__init__(self, data=alongpath[:, 1:],
                                 columns=header[1:], geometry=pts)

        self.path = LineString(pts)

        self['DateTime'] = datetime
        self.set_index('Timestep', inplace=True, drop=False)

    def calculate_vector(self, reverse=False):
        """
        Calculate the following in radians:
            -The bearings from origin to each timestep
            -The bearings between timesteps (closer to farther from origin)
            -The circular mean of the origin-timestep bearings

        Each timestep contains the bearing needed to get to that timestep from
        the origin (``bearings_from_origin`` or ``bearings_from_origin_r``)
        or to get to that timestep from the previous timestep (``bearings_ptp``
        or ``bearings_ptp_r``)

        Parameters
        ----------
        reverse : Boolean
            Default False.  Indicates the original trajectory or its reversed
            counterpart.  Reversed trajectory must first be loaded via
            ``self.load_reversetraj()``

        """
        labels = {False: ['bearings_from_origin', 'bearings_ptp',
                          'circular_mean'],
                  True: ['bearings_from_origin_r', 'bearings_ptp_r',
                         'circular_mean_r']}

        which_traj = {False: 'path',
                      True: 'path_r'}

        try:
            lat, lon = np.radians(getattr(self, which_traj[reverse]).xy)
        except:
            raise AttributeError('Reversed trajectory is not loaded!')

        a = np.cos(lat) * np.sin(lon - lon[0])
        b = (math.cos(lat[0]) * np.sin(lat) -
             math.sin(lat[0]) * np.cos(lat) * np.cos(lon - lon[0]))

        bearings_fo = np.arctan2(a, b)

        # Set bearings from origin column
        self[labels[reverse][0]] = bearings_fo

        x = np.mean(np.cos(bearings_fo))
        y = np.mean(np.sin(bearings_fo))

        # Calculate circular means
        setattr(self, labels[reverse][2], math.atan2(y, x))

        a = np.cos(lat[1:]) * np.sin(lon[1:] - lon[:-1])
        b = (np.cos(lat[:-1]) * np.sin(lat[1:]) -
             np.sin(lat[:-1]) * np.cos(lat[1:]) * np.cos(lon[1:] - lon[:-1]))

        bearings_ptp = np.arctan2(a, b)

        # point to point bearings column, first entry is 0
        self[labels[reverse][1]] = 0.0
        self.loc[self.index[1:], labels[reverse][1]] = bearings_ptp

    def calculate_distance(self, reverse=False):
        """
        Calculate the following great circle distances in meters:
            -Between each timepoint
            -Cumulative along-path travel distance from origin
            -Distance between origin and point

        Parameters
        ----------
        reverse : Boolean
            Default False.  Indicates the original trajectory or its reversed
            counterpart.  Reversed trajectory must first be loaded via
            ``self.load_reversetraj()``.

        """

        which_traj = {False: 'path',
                      True: 'path_r'}

        labels = {False: ['Distance_ptp', 'Cumulative_Dist',
                          'Dist_from_origin'],
                  True: ['Distance_ptp_r', 'Cumulative_Dist_r',
                         'Dist_from_origin_r']}

        lat, lon = np.radians(getattr(self, which_traj[reverse]).xy)

        distance = np.empty((lat.size))

        distance[0] = 0.0
        distance[1:] = (np.arccos(np.sin(lat[1:]) * np.sin(lat[:-1]) +
                                  np.cos(lat[1:]) * np.cos(lat[:-1]) *
                                  np.cos(lon[:-1] - lon[1:])) * 6371) * 1000

        self[labels[reverse][0]] = distance

        self[labels[reverse][1]] = np.cumsum(distance)

        dist = (np.arccos(np.sin(lat[1:]) * np.sin(lat[0]) +
                          np.cos(lat[1:]) * np.cos(lat[0]) *
                          np.cos(lon[0] - lon[1:])) * 6371) * 1000

        self[labels[reverse][2]] = 0.0
        self.loc[self.index[1:], labels[reverse][2]] = dist

    def distance_between2pts(self, coord0, coord1):
        """
        Calculate distance between two sets of coordinates

        Parameters
        ----------
        coord0 : tuple of floats
            Latitude, longitude coordinate pair in degrees
        coord1 : tuple of floats
            Latitude, longitude coordinate pair in degrees

        Returns
        -------
        Distance : float
            Great circle distance in meters.

        """
        coord0 = np.radians(coord0)
        coord1 = np.radians(coord1)

        distance = (np.arccos(np.sin(coord1[0]) * np.sin(coord0[0]) +
                              np.cos(coord1[0]) * np.cos(coord0[0]) *
                              np.cos(coord0[1] - coord1[1])) * 6371) * 1000

        return distance

    def find_destination(self, lat0, lon0, bearing, distance):
        """
        Find the destination given a bearing and a set of coordinates.

        Parameters
        ----------
        lat0 : float
            Latitude of starting point in degrees
        lon0 : float
            Longitude of starting point in degrees
        bearing : float
            Direction from starting point in radians (``self.circular_mean``
            is in degrees)

        Returns
        -------
        latx : float
            Latitude of destination in degrees
        lonx : float
            Longitude of destination in degrees

        """

        d2r = distance / 6371000

        latr = math.radians(lat0)
        lonr = math.radians(lon0)

        latx = math.asin(math.sin(latr) * math.cos(d2r) +
                         math.cos(latr) * math.sin(d2r) *
                         math.cos(bearing))

        lonx = math.degrees(lonr + math.atan2(math.sin(bearing) *
                                              math.sin(d2r) * math.cos(latr),
                                              math.cos(d2r) - math.sin(latr) *
                                              math.sin(latx)))
        latx = math.degrees(latx)

        return latx, lonx
