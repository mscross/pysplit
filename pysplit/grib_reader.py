from __future__ import division
import numpy as np
import gdal
import gdalconst
import fnmatch
import pyhdf.SD
import itertools
import os
import mapmaker as mm


def get_bandnum(level, variable, datafile):
    """
    Reads GRIB (ncar) files, returns number of desired raster band.

    Parameters
    ----------
    level : string
        Indicates which pressure level or surface the reader should search for
    variable : string
        Indicates the variable that the reader is searching for
    datafile : string
        The file that the reader should look through

    Returns
    -------
    bandnum : int
        Number of the raster band required

    """

    # Open the datafile
    dataset = gdal.Open(datafile, gdalconst.GA_ReadOnly)

    bandnum = 0

    for i in range(1, dataset.RasterCount + 1):
        # Open band and get metadata
        band = dataset.GetRasterBand(i)
        metadata = band.GetMetadata()

        # Get the band surface level and the variable
        bandlevel = metadata['GRIB_SHORT_NAME']
        bandvar = metadata['GRIB_ELEMENT']

        # Break if this level and variable is the one we're looking for
        if (variable == bandvar) and (level == bandlevel):
            bandnum = i
            break

    if bandnum == 0:
        raise ValueError('Band not found!  Check surface, band variable using',
                         ' getsurface_getvar().')

    return bandnum


def getsurface_getvar(datafile):
    """
    Opens GRIB (ncar) files, prints the surface level
        and variable of each raster band.

    Parameters
    ----------
    datafile : string
        The file whose metadata is retrieved

    """

    dataset = gdal.Open(datafile, gdalconst.GA_ReadOnly)

    if dataset.RasterCount == 444 :
        stop_here = 37 + 1
    else:
        stop_here = dataset.RasterCount + 1

    print 'GRIB_SHORT_NAME \t GRIB_ELEMENT'

    for i in range(1, stop_here):
        band = dataset.GetRasterBand(i)
        metadata = band.GetMetadata()

        # Get the band surface level and the variable
        bandlevel = metadata['GRIB_SHORT_NAME']
        bandvar = metadata['GRIB_ELEMENT']

        print bandlevel + '\t\t\t\t' + bandvar

    return None


def file_lister(filedates, startstring, data_dir):
    """
    Get a list of files that have both 1 of a series of dates
        and match startstring.

    Parameters
    ----------
    filedates : list of strings
        List of 6 character strings representing years and months
        of interest
    startstring : string
        Initial part of file names
    data_dir : string
        Full or relative path to the data directory

    Returns
    -------
    completefilenames : list of files
        List of files that have both startstring and one of the
        dates of interest

    """

    # Initialize
    orig_dir = os.getcwd()
    filelist = []

    try:
        os.chdir(data_dir)

        # Get list of files
        _, _, files = next(os.walk('.'))

        # Find matching files
        for eachfile, date in itertools.product(files, filedates):
            if fnmatch.fnmatch(eachfile, startstring + '*' + date + '*'):
                filelist.append(eachfile)

    finally:
        os.chdir(orig_dir)

    completefilenames = []

    for f in filelist:
        completefilenames.append(os.path.join(data_dir, f))

    return completefilenames


def file_dates(years, years_isrange, months, months_isrange):
    """
    Constructs list of YYYYMM date strings from the given
        years and months of interest

    Parameters
    ----------
    years : list of ints
        The years of interest
    years_isrange : Boolean
        Indicates whether to use `years` as is or to construct
        a range of years between years[0] and years[-1]
    months : list of ints
        The months of interest
    months_isrange : Boolean
        Indicates whether to use `months` as is or to construct
        a range of months between months[0] and months[-1]

    Returns
    -------
    dates : list of strings
        List of strings representing dates in the format YYYYMM

    """

    dates = []

    if years_isrange:
        years = range(years[0], years[-1] + 1)
    if months_isrange:
        months = range(months[0], months[-1] + 1)

    for y, m in itertools.product(years, months):
        dates.append('{0:04}{1:02}'.format(y, m))

    return dates


def get_gribdata(level, var, startstring, data_dir, filedates,
                 filelist=None):
    """
    Gets a list of GRIB files, determines which band contains the desired data.
        For each file, the data is read into an array, and the
        data array is put into a list.  Also returns latitudes and longitudes
        that represent the data grid.

    Parameters
    ----------
    level : string
        Desired surface level of data
    var : string
        Band variable
    startstring : string
        Initial part of file names.  Only needed if filelist is None.
    data_dir : string
        Full or relative path to the data directory.  Only needed if
        filelist is None.
    filedates : list of strings
        Default None. List of strings representing dates in the format YYYYMM.
        Generated from file_dates().  Used to gather a list of filenames.
        Only needed if filelist is None.
    filelist : list of strings
        A list of filenames.

    Returns
    -------
    gribdata : list of (M, N) ndarray of floats
        Arrays of data from all GRIB files in indicated time period
    latitudes : ndarray of floats
        (M) array of latitudes of gribdata
    longitudes : ndarray of floats
        (N) array of longitudes of gribdata

    """

    # Create filelist from startstring and filedates if not given
    if filelist is None:
        filelist = file_lister(filedates, startstring, data_dir)

    gribdatalist = []

    # Get the number of the band that contains the desired data
    bandnum = get_bandnum(level, var, filelist[0])

    # Open each file, get band, read band
    # into an array and put array into a list
    for files in filelist:

        dataset = gdal.Open(files, gdalconst.GA_ReadOnly)

        band = dataset.GetRasterBand(bandnum)

        bandval = band.ReadAsArray(0, 0, band.XSize, band.YSize)

        gribdatalist.append(bandval)

    # Get the spatial coverage and start/end of the data
    geotransform = dataset.GetGeoTransform(0)

    # Generate longitudes and latitudes
    longitudes = np.linspace(geotransform[0], 360.0 + geotransform[0],
                             band.XSize, endpoint=False)
    latitudes = np.linspace(geotransform[3], -1.0 * geotransform[3],
                            band.YSize)

    return gribdatalist, latitudes, longitudes


def get_hdfdata(filedates, startstring, data_dir, var='precipitation',
                filelist=None):
    """
    Gets a list of HDF files and for each file, reads the data in an array.
        The data arrays are put into a list.  Also returns the latitudes and
        longitudes that represent the data grid

    Parameters
    ----------
    filedates : list of strings
        Default None. List of strings representing dates in the format YYYYMM.
        Generated from file_dates().  Used to gather a list of files.
    startstring : string
        Initial part of file names
    data_dir : string
        Path to data location
    var : string
        Default 'precipitation'.  The variable to select from the HDF file
    filelist : list of strings
        Default None.  The list of HDF files to introspect.  If ``None``, list
        will be assembled from ``filedates`` and ``startstring``, otherwise
        ``filelist`` will override args.

    Returns
    -------
    hdfdata : list of (M, N) ndarray of floats
        arrays of data from all trmm (HDF) files given
    longitudes : (M) ndarray of floats
        array of longitudes of hdfdata
    latitudes : (N) ndarray of floats
        array of latitudes of hdfdata

    """
    if filelist is None:
        filelist = file_lister(filedates, startstring, data_dir)

    hdfdatalist = []

    for files in filelist:

        hdf_file = pyhdf.SD.SD(files)

        variable = hdf_file.select(var)

        vardata = variable.get()

        hdfdatalist.append(vardata)

    # Set for trmm files
    longitudes = np.arange(-180.0, 180.0, 0.25)
    latitudes = np.arange(-49.875, 50.125, 0.25)

    return hdfdatalist, longitudes, latitudes


def averager(gridded_data):
    """
    Averages list of numpy arrays along 2nd axis

    Parameters
    ----------
    gridded_data : list of (M, N) ndarrays of floats

    Returns
    -------
    gridded_data : (M, N) ndarray of floats

    """

    gridded_data = np.mean(np.dstack(gridded_data), axis=2)

    return gridded_data


def windcontours(windmap, uband, vband, latitudes, longitudes, contourf=True,
                 limits=None, **kwargs):
    """
    Contour map of windspeed.

    Parameters
    ----------
    windmap : Basemap instance
        Initialize a basemap first using MapDesign.make_basemap()
    uband : (M, N) ndarray of floats
        Array of the U component of wind velocity
    vband : (M, N) ndarray of floats
        Array of the V component of wind velocity
    latitudes : (M) ndarray of floats
        Array of latitudes of uband, vband
    longitudes : (N) ndarray of floats
        Array of longitudes of uband, vband

    Keyword Arguments
    -----------------
    limits : list of floats or ints
        Default None.  Used to trim zdata grid, longitudes, latitudes.
        [bottom latitude, left longitude, top latitude, right longitude]

    Other Parameters
    ----------------
    kwargs : passed to mm.meteo_contouring(), then Basemap.contour() or
        Basemap.contourf(), then the axis methods

    Returns
    -------
    w : matplotlib PathCollection instance
        Mappable for use in creating colorbars.  Colorbars may be created
        using make_cbar() or make_cax_cbar().

    """

    # Trim grid to a more localized box
    if limits is not None:
        _, _, uband = gridlimit(limits, longitudes, latitudes, uband)
        longitudes, latitudes, vband = gridlimit(limits, longitudes, latitudes)

    # Calculate windspeed
    windspeed = np.sqrt(uband ** 2 + vband ** 2)

    # Plot data as filled contours
    w = mm.meteo_contouring(windmap, windspeed, longitudes, latitudes,
                            contourf=contourf, **kwargs)

    return w


def plot_othervar(varmap, vardata, latitudes, longitudes, contourf=True,
                  is_geopotential=False, filetype='ncar', limits=None,
                  **kwargs):
    """
    For plotting any variable that is not wind.  Interrogates ncar (grib) or
        trmm (hdf) files on any surface and plots data from specified years or
        interval of years in a filled contour plot on a basemap.

    Parameters
    ----------
    varmap : Basemap instance
        Initialize a basemap first using MapDesign.make_basemap()
    vardata : (M, N) ndarray of floats
        Data array
    latitudes : (M) ndarray of floats
        Array of latitudes of vardata
    longitudes : (N) ndarray of floats
        Array of longitudes of vardata

    Keyword Arguments
    -----------------
    filetype : string
        Default 'ncar'.  ['ncar'|'trmm'].  TRMM data has a different
        orientation than NCAR data and must be transposed before plotting.
    limits : list of floats or ints
        Default None.  Used to trim zdata grid, longitudes, latitudes.
        [bottom latitude, left longitude, top latitude, right longitude]

    Other Parameters
    ----------------
    kwargs : passed to mm.meteo_contouring(), then Basemap.contour() or
        Basemap.contourf(), then the axis methods

    Returns
    -------
    v : matplotlib PathCollection instance
        Mappable for use in creating colorbars.  Colorbars may be created
        using make_cbar() or make_cax_cbar().

    """

    # TRMM files are oriented weirdly
    if filetype is 'trmm':
        vardata = vardata.T

    # Trim grid to a more localized box
    if limits is not None:
        longitudes, latitudes, vardata = gridlimit(limits, longitudes,
                                                   latitudes, vardata)

    # Calculate meters from geopotential
    if is_geopotential:
        vardata = vardata / 9.80665

    # Make the contour plot
    v = mm.meteo_contouring(varmap, vardata, longitudes, latitudes,
                            contourf=contourf, **kwargs)

    return v


def windbarbs(windmap, uband, vband, latitudes, longitudes, vectors='arrows',
              limits=None, zorder=20, scale=1000, **kwargs):
    """
    Map wind velocity using vectors (arrows or barbs).

    plt.quiverkey(q, 0.25, 0.05, 60, '60 m/s', labelpos='W',
              labelcolor='red', fontproperties={'size': 20},
              zorder=zorder, color='red')

    Parameters
    ----------
    windmap : Basemap instance
        Initialize a basemap first using MapDesign.make_basemap()
    uband : (M, N) ndarray of floats
        Array of the U component of wind velocity
    vband : (M, N) ndarray of floats
        Array of the V component of wind velocity
    latitudes : (M) ndarray of floats
        Array of latitudes of uband, vband
    longitudes : (N) ndarray of floats
        Array of longitudes of uband, vband
    vectors : string
        Default 'arrows'.  ['arrows'|'barbs']
        The vector type to indicate wind velocity.  'arrows' comes with a key.
    limits : list of floats or ints
        Default None.  Used to trim zdata grid, longitudes, latitudes.
        [bottom latitude, left longitude, top latitude, right longitude]
    zorder : int
        Default 20.  The zorder of the vectors `windmap`
    scale : int
        Default 1000. Used to adjust arrow size.

    Other Parameters
    ----------------
    kwargs : passed to Basemap.quiver() or Basemap.barbs(), then the axis
        methods

    Returns
    -------
    q : collection
        Collection of barbs or arrows

    """

    # Trim grid to a more localized box
    if limits is not None:
        _, _, uband = gridlimit(limits, longitudes, latitudes, uband)
        longitudes, latitudes, vband = gridlimit(limits, longitudes, latitudes)

    # Change longitudes to -180 to 180
    longitudes = np.where(longitudes > 180.0, longitudes - 360.0, longitudes)

    longitude_sort = np.argsort(longitudes)

    longitudes = longitudes[longitude_sort]
    uband = uband[:, longitude_sort]
    vband = vband[:, longitude_sort]

    # Transform vectors to lie correctly on projection
    u, v, lns, lts = windmap.transform_vector(uband[::-1, :], vband[::-1, :],
                                              longitudes, latitudes[::-1],
                                              12, 16, returnxy=True)
    # Plot vectors as arrows or barbs
    if vectors is 'arrows':
        q = windmap.quiver(lns, lts, u, v, zorder=zorder, scale=1000, **kwargs)

    elif vectors is 'barbs':
        q = windmap.barbs(lns, lts, u, v, zorder=zorder, **kwargs)

    return q


def gridlimit(limits, longitudes, latitudes, data):
    """
    Get rid of portion of grid beyond the given limits

    Parameters
    ----------
    limits : list of floats or ints
        [bottom latitude, left longitude, top latitude, right longitude]
    longitudes : (N) ndarray of floats
        Array of longitudes of data
    latitudes : (M) ndarray of floats
        Array of latitudes of data
    data : (M, N) ndarray of floats
        Data array

    Returns
    -------
    longitudes : (S) ndarray of floats
        Potentially trimmed array of longitudes of data, now in -180 to 180.
    latitudes : (R) ndarray of floats
        Potentially trimmed array of latitudes of data
    data : (R, S) ndarray of floats
        Trimmed data array

    """

    if limits[0] is not None:
        lat_inds = np.nonzero(latitudes > limits[0])[0]
        latitudes = latitudes[lat_inds]
        data = data[lat_inds, :]

    if limits[2] is not None:
        lat_inds = np.nonzero(latitudes < limits[2])[0]
        latitudes = latitudes[lat_inds]
        data = data[lat_inds, :]

    if limits[1] > limits[3] and limits[1] > 180:
        lon_tmp = []
        for lon in longitudes:
            if lon > 180:
                lon = 360 - lon
            lon_tmp.append(lon)
        longitudes = np.asarray(lon_tmp)

    if limits[1] is not None:
        lon_inds = np.nonzero(longitudes > limits[1])[0]
        longitudes = longitudes[lon_inds]
        data = data[:, lon_inds]

    if limits[3] is not None:
        lon_inds = np.nonzero(longitudes < limits[3])[0]
        longitudes = longitudes[lon_inds]
        data = data[:, lon_inds]

    return longitudes, latitudes, data
