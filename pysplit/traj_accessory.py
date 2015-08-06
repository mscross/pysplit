from __future__ import division
import numpy as np
import math


def tracemean_vector(latitude, longitude):
    """
    Calculates the bearing of the mean vector of a trace with points
    (``latitude``, ``longitude``)

    Parameters
    ----------
    latitude : iterable of floats
        The latitudes of a set of points
    longitude ; iterable of floats
        The longitudes of a set of points

    Returns
    -------
    mean_vector : float
        The mean bearing of the trace
    bearings : 1D ndarray of floats
        The bearing at t0 along great circles between t0 and
        each point of the trace

    """

    bearings = []

    # Initialize the first point of the set of points
    t0 = (latitude[0], longitude[0])

    # For each point, find the bearing at t0 along a great circle
    # between t0 and tx
    for lat, lon in zip(latitude, longitude):
        tx = (lat, lon)
        bearings.append(great_circle_bearing(t0, tx))

    bearings = np.asarray(bearings).astype(np.float64)

    # Find the circular mean of the bearings
    mean_vector = circular_means(bearings)

    return mean_vector, bearings


def great_circle_bearing(t0, tx):
    """
    Calculates bearing at point t0 along great circle towards point tx

    Parameters
    ----------
    t0 : tuple of floats
        (Latitude, longitude) coordinate pair.  Bearing calculated at
        this point (bearing along great circle changes with position)
    tx : tuple of floats
        (Latitude, longitude) coordinate pair; endpoint of great circle.

    Returns
    -------
    bearing : float
        Bearing in degrees at t0 along a great circle between points t0, tx

    """

    # Convert degrees to radians
    lat0 = math.radians(t0[0])
    lon0 = math.radians(t0[1])
    latx = math.radians(tx[0])
    lonx = math.radians(tx[1])

    a = math.cos(latx) * math.sin(lonx - lon0)
    b = (math.cos(lat0) * math.sin(latx) -
         math.sin(lat0) * math.cos(latx) * math.cos(lonx - lon0))

    bearing = math.degrees(math.atan2(a, b))

    return bearing


def circular_means(bearings):
    """
    Average an array of bearings by finding the circular mean

    The circular mean is the correct way to average angles.  For example,
    the arithmetic mean of a vector at 90 degrees and one at 270 (-90)
    degrees is 180 (0) degrees (i.e., the average of east and west is south
    or north).  The circular mean correctly calculates that the average is
    a vector of length 0 with no direction.

    Parameters
    ----------
    bearings : 1D ndarray of floats
        Array of (great circle) bearings from the location of time 0 to the
        location of each point along a trace

    Returns
    -------
    circ_mean : float
        The average bearing, found by taking the circular mean

    """

    x = []
    y = []

    # Convert from polar coordinates on unit circle to Cartesian coordinates
    for b in bearings:
        b = math.radians(b)
        x.append(math.cos(b))
        y.append(math.sin(b))

    x = np.asarray(x).astype(np.float64)
    y = np.asarray(y).astype(np.float64)

    # Take arithmetic mean of Cartesian coordinates
    x_mean = np.mean(x)
    y_mean = np.mean(y)

    # Convert back to polar coordinates
    # Point will be on unit disk; don't care about r
    circ_mean = math.degrees(math.atan2(y_mean, x_mean))

    return circ_mean


def distance_overearth(latitude, longitude):
    """
    Calculate the distance between points on the earth's surface

    Parameters
    ----------
    latitude : iterable of floats of at least length 2
        The latitudes of a set of points on the surface of the earth
    longitude : iterable of floats of at least length 2
        The longitudes of a set of points on the surface of the earth

    Returns
    -------
    distance : 1D ndarray of floats
        The distance in meters between each latitude, longitude point
        and the previous latitude, longitude point.  The first item in
        distance is 0, there being no previous point.

    """

    # Convert to radians
    lat_rad = np.radians(latitude)
    lon_rad = np.radians(longitude)

    # Offset latitude and longitude each by 1
    lat0 = lat_rad[:-1]
    lat1 = lat_rad[1:]

    lon0 = lon_rad[:-1]
    lon1 = lon_rad[1:]

    distance = (np.arccos(np.sin(lat1) * np.sin(lat0) +
                          np.cos(lat1) * np.cos(lat0) *
                          np.cos(lon0 - lon1)) * 6371) * 1000

    # Pad array so first element of distance is 0
    distance = np.pad(distance, (1, 0), 'constant',
                      constant_values=(0.0, 0.0))

    return distance


def sum_distance(distance):
    """
    Take a series of distances between points along a trace and find the
        cumulative distance between each point and time 0.

    Parameters
    ----------
    distance : 1D ndarray of floats
        The distance in meters between a point and a previous point
        along a trace on the earth's surface.

    Returns
    -------
    total_distance : 1D ndarray of floats
        The distance between a point in a trace on the earth's surface
        and the point at time 0 of the trace.

    """

    tot_dist = []

    # At each point, add up the distance between each prior step
    for i in range(0, distance.size):
        dist = sum(distance[: i + 1])
        tot_dist.append(dist)

    total_distance = np.asarray(tot_dist).astype(np.float64)

    return total_distance


def find_destination(t0, bearing, distance, unit='m'):
    """
    Calculate destination given initial location, bearing, and distance.

    Bearing is taken at t0 along great circle arc from t0 to destination.
    Related to ``distance_overearth()``, ``great_circle_bearing()``.

    Parameters
    ----------
    t0 : tuple of floats
        (Latitude, Longitude) coordinates of starting point
    bearing : float
        The bearing in degrees at ``t0`` along a great circle between
        t0 and the destination
    distance : float
        The distance in meters, kilometers, or miles
        between ``t0`` and the destination
    unit : string
        Default 'm'.  ['m'|'km'|'mi'].  The unit of ``distance``

    Returns
    -------
    tx : tuple of floats
        (Latitude, Longitude) coordinates of destination at ``bearing``
        and ``distance`` relative to ``t0``

    """

    meters = ['m', 'meter', 'meters']
    kilometers = ['km', 'kilometer', 'kilometers']

    # Initialize earth radius in appropriate units
    if unit in meters:
        earth_r = 6371 * 1000
    elif unit in kilometers:
        earth_r = 6371
    else:
        earth_r = 6371
        distance = convert_mi2km(distance)
        print 'Converted distance in miles to kilometers: ' + distance

    bearing = math.radians(bearing)
    lat0 = t0[0]
    lon0 = t0[1]

    # Calculate new lat, lon
    latx = math.asin(math.sin(lat0) * math.cos(distance / earth_r) +
                     math.cos(lat0) * math.sin(distance / earth_r) *
                     math.cos(bearing))

    lonx = lon0 + math.atan2(math.sin(bearing) * math.sin(distance / earth_r) *
                             math.cos(lat0),
                             math.cos(distance / earth_r) -
                             math.sin(lat0) * math.sin(latx))

    tx = (latx, lonx)

    return tx


def convert_mi2km(distance):
    """
    Convert miles to kilometers.

    Parameters
    ----------
    distance : float or iterable of floats
        A value or values in units of miles

    Returns
    -------
    distance_km : float
        The given value or values converted to kilometers

    """

    distance_km = distance / 0.62137

    return distance_km


def convert_km2mi(distance):
    """
    Convert kilometers to miles.

    Parameters
    ----------
    distance : float or iterable of floats
        A value or values in units of kilometers

    Returns
    -------
    distance_mi : float
        The given value or values converted to miles

    """

    distance_mi = distance * 0.62137

    return distance_mi


def convert_w2q(mixing_ratio):
    """
    Convert mixing ratio to specific humidity.

    Parameters
    ----------
    mixing_ratio : iterable of floats, or float
        g water per kg of air

    Returns
    -------
    q : iterable of floats, or float
        Specific humidity/humidities
        g water per kg of dry air

    """

    # Convert from g/kg to kg/kg to g/kg
    mr_kg = mixing_ratio / 1000
    q_kg = mr_kg / (mr_kg + 1)
    q = q_kg * 1000

    return q


def convert_q2w(specific_humidity):
    """
    Convert specific humidity to mixing ratio

    Parameters
    ----------
    specific_humidity : iterable of floats, or float
        g or kg water per kg of dry air

    Returns
    -------
    mixing_ratio : iterable of floats, or float
        g or kg water per kg of air

    """

    # Convert from g/kg to kg/kg to g/kg
    sph_kg = specific_humidity / 1000
    w_kg = sph_kg / (1 - sph_kg)
    w = w_kg * 1000

    return w


def convert_w2rh(mixing_ratio, temperature, pressure):
    """
    Convert ``mixing ratio`` (absolute humidity) to relative humidity

    Parameters
    ----------
    mixing_ratio : iterable of floats, or float
        g or kg water per kg of air
    temperature : iterable of floats, or float
        Temperature in K or degrees C.  Automatic conversion of units,
        for temperatures reasonable at the earth's surface
    pressure : iterable of floats, or float
        Atmospheric pressure in hectopascals

    Returns
    -------
    rh : iterable of floats, or float
        relative humidity on a scale of 0.0 to 100.0 %

    """

    try:
        relhumid_ls = []

        for w, t, p in zip(mixing_ratio, temperature, pressure):
            # Detect if temperature is in K or degrees C, convert to degrees C
            if t > 100.0:
                t_celsius = t - 273.15
            else:
                t_celsius = t

            # Calculate saturation vapor pressure, saturation mixing ratio
            satvapor = 6.11 * (10.0 ** ((7.5 * t_celsius) /
                                        (237.7 + t_celsius)))
            sat_w = 621.97 * (satvapor / (p - satvapor))

            rh = (w / sat_w) * 100.0
            relhumid_ls.append(rh)

        rh = np.asarray(relhumid_ls).astype(np.float64)

    except:
        # Detect if temperature is in K or degrees C, convert to degrees C
        if temperature > 100.0:
            t_celsius = temperature - 273.15
        else:
            t_celsius = temperature

        # Calculate saturation vapor pressure, saturation mixing ratio
        satvapor = 6.11 * (10.0 ** ((7.5 * t_celsius) / (237.7 + t_celsius)))
        sat_w = 621.97 * (satvapor / (pressure - satvapor))

        rh = (mixing_ratio / sat_w) * 100.0

    return rh


def convert_rh2w(relative_humidity, temperature, pressure):
    """
    Convert ``relative_humidity`` to mixing ratio (absolute humidity)

    Parameters
    ----------
    relative_humidity : iterable of floats, or float
        Relative humidity on a scale of 0.0 to 100.0 %
    temperature : iterable of floats, or float
        Temperature in K or degrees C.  Automatic conversion of units,
        for temperatures reasonable at the earth's surface
    pressure : iterable of floats, or float
        Atmospheric pressure in hectopascals

    Returns
    -------
    mixing_ratio : iterable of floats, or float
        g or kg water per kg of air

    """

    try:
        w_ls = []

        for rh, t, p in zip(relative_humidity, temperature, pressure):
            # Detect if temperature is in K or degrees C, convert to degrees C
            if t > 100.0:
                t_celsius = t - 273.15
            else:
                t_celsius = t

            # Calculate saturation vapor pressure, saturation mixing ratio
            satvapor = 6.11 * (10.0 ** ((7.5 * t_celsius) /
                                        (237.7 + t_celsius)))
            sat_mixratio = 621.97 * (satvapor / (p - satvapor))

            w = (rh / 100.0) * sat_mixratio
            w_ls.append(w)

        w = np.asarray(w_ls).astype(np.float64)

    except:
        # Detect if temperature is in K or degrees C, convert to degrees C
        if temperature > 100.0:
            t_celsius = temperature - 273.15
        else:
            t_celsius = temperature

        # Calculate saturation vapor pressure, saturation mixing ratio
        satvapor = 6.11 * (10.0 ** ((7.5 * t_celsius) / (237.7 + t_celsius)))
        sat_mixratio = 621.97 * (satvapor / (pressure - satvapor))

        w = (relative_humidity / 100.0) * sat_mixratio

    return w


def geographic_midpt(t0, tx):
    """
    Calculate the midpoint between two coordinate pairs on the Earth's surface.

    The geographic midpoint is calculated by finding the center of gravity
        between two locations.

    Parameters
    ----------
    t0 : tuple of floats
        Latitude, longitude (in decimal degrees) of a point
    tx : tuple of floats
        Latitude, longitude (in decimal degrees) of a second point

    Returns
    -------
    lat : float
        The latitude of the midpoint
    lon : float
        The longitude of the midpoint

    """

    # Convert from degrees to radians
    lat1 = t0[0] * (math.pi / 180)
    lon1 = t0[1] * (math.pi / 180)
    lat2 = tx[0] * (math.pi / 180)
    lon2 = tx[1] * (math.pi / 180)

    # Convert lat/lon to Cartesian coordinates
    x1 = math.cos(lat1) * math.cos(lon1)
    y1 = math.cos(lat1) * math.sin(lon1)
    z1 = math.sin(lat1)

    x2 = math.cos(lat2) * math.cos(lon2)
    y2 = math.cos(lat2) * math.sin(lon2)
    z2 = math.sin(lat2)

    # Find the average
    x = (x1 + x2) / 2
    y = (y1 + y2) / 2
    z = (z1 + z2) / 2

    # Find average longitude and latitude in radians
    lon = math.atan2(y, x)
    hyp = math.sqrt(x * x + y * y)
    lat = math.atan2(z, hyp)

    # Convert to degrees
    lat = lat * (180 / math.pi)
    lon = lon * (180 / math.pi)

    return lat, lon


def grid_data(x, y, data, cell_value, binsize):
    """
    Place unevenly spaced 2D data on a grid by 2D binning using nearest
    neighbor interpolation.

    Originally Example 3 by ccampo (2010-07-11) from
    http://wiki.scipy.org/Cookbook/Matplotlib/Gridding_irregularly_spaced_data,

    Adapted for PySPLIT by mscross, 2014-06-12.

    Parameters
    ----------
    x : 1D ndarray of scalars
        The independent data x-axis of the grid (longitude)
    y : 1D ndarray of scalars
        The independent data y-axis of the grid (latitude)
    z : 1D ndarray of scalars
        The unevenly spaced dependent data.
        For example, specific humidity at each x,y point along a trajectory
    cell_value : string
        Determines the value of each cell from the contents of the bin.
        ['median'|'mean'|'cumulative'|'max'|'min'|'range'|'stdev']
    binsize : float
        The width, height of each bin.  Only square bins supported.

    Returns
    -------
    grid : masked 2D ndarray of scalars
        The evenly gridded data.  The value of each cell in relation to the
        bin contents is determined by ``cell_value``.
        Invalid values are masked.
    xi : 2D ndarray of floats
        The grid of x bin bounds
    yi : 2D ndarray of floats
        The grid of y bin bounds
    bins : 2D ndarray of floats
        A grid the same shape as ``grid``, except the value of each cell is
        the number of points in the bin.
    wherebin : 2D list
        A 2D list the same shape as ``grid`` and ``bins`` where each cell
        contains the indices of ``data`` that correspond to the values stored
        in the particular bin.

    """

    # Initialize dictionary
    cell_value_dict = {'cumulative' : np.sum,
                       'mean' : np.mean,
                       'median' : np.median,
                       'max' : np.max,
                       'min' : np.min,
                       'stdev' : np.std,
                       'range' : np.ptp}

    # Get extreme longitudes and latitudes
    xmin = x.min()
    xmax = x.max()

    ymin = y.min()
    ymax = y.max()

    # Make coordinate arrays
    xi = np.arange(xmin, xmax + binsize, binsize)
    yi = np.arange(ymin, ymax + binsize, binsize)

    xi, yi = np.meshgrid(xi, yi)

    # Make `grid`
    grid = np.zeros(xi.shape, dtype=x.dtype)
    nrow, ncol = grid.shape

    # Set up `bins` and `wherebin` if either/both are called for
    bins = np.copy(grid)
    wherebin = np.copy(grid)
    wherebin = wherebin.tolist()

    # Fill in the grid
    for row in range(nrow):
        for col in range(ncol):

            # Get x, y coordinates of current position
            xc = xi[row][col]
            yc = yi[row][col]

            # Find the position(s) in the original data array that
            # xc, yc correspond to
            # Get absolute values of all items in x, y - xc, yc
            posx = np.abs(x - xc)
            posy = np.abs(y - yc)

            # Get boolean array of where condition is met in both x, y
            # ibin.size == x.size == y.size
            ibin = np.logical_and(posx < binsize / 2.0, posy < binsize / 2.0)

            # Get array of indices of data points that are in the current pos
            ind = np.where(ibin)[0]

            # Fill bin (True values in ibin will put data value into bin)
            bin = data[ibin]

            # Set cell values in grid
            if bin.size != 0:
                binval = cell_value_dict[cell_value](bin)
                grid[row, col] = binval
            else:
                grid[row, col] = -999.0

            # Update `wherebin` and `bins` if necessary
            wherebin[row][col] = ind
            bins[row, col] = bin.size

    # Mask 'invalid' entries.  PySPLIT follows the convention of
    # -999.0 as a fill value for invalid or missing data
    grid = np.ma.masked_less_equal(grid, -999.0)

    return grid, xi, yi, bins, wherebin
