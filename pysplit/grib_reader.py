import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import shiftgrid
import gdal
import ogr
from gdalconst import *
import mapmaker as mm
import fnmatch
import pyhdf.SD
import itertools
import os


def getband_latslons(level, variable, datafile):
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
    latitudes : ndarray of floats
        The latitudes of the data grid
    longitudes : ndarray of floats
        The longitudes of the data grid

    """

    # Open the datafile
    dataset = gdal.Open(datafile, GA_ReadOnly)

    bandnum = 0

    for i in range (1, dataset.RasterCount + 1):
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

    # Get the spatial coverage and start/end of the data
    geotransform = dataset.GetGeoTransform(0)

    # Generate longitudes and latitudes
    longitudes = np.linspace(geotransform[0], 360.0 + geotransform[0],
                             band.XSize, endpoint=False)
    latitudes = np.linspace(geotransform[3], -1.0 * geotransform[3], band.YSize)

    return bandnum, latitudes, longitudes



def getsurface_getvar(datafile):
    """
    Opens GRIB (ncar) files, prints the surface level
        and variable of each raster band.

    Parameters
    ----------
    datafile : string
        The file whose metadata is retrieved

    """

    dataset = gdal.Open(datafile, GA_ReadOnly)

    if dataset.RasterCount == 444 :
        stop_here = 37 + 1
    else:
        stop_here = dataset.RasterCount + 1

    print 'GRIB_SHORT_NAME \t GRIB_ELEMENT'

    for i in range (1, stop_here):
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
            matchstring = startstring + '*' + date + '*'
            if fnmatch.fnmatch(eachfile, startstring +'*'+date+'*'):
                filelist.append(eachfile)

    finally:
        os.chdir(orig_dir)

    completefilenames=[]

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
        years = range(years[0], years[-1]+1)
    if months_isrange:
        months = range(months[0], months[-1]+1)

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

    Keyword Arguments
    -----------------
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
    bandnum, latitudes, longitudes = getband_latslons(level, var, filelist[0])

    # If the appropriate band was found, open each file, get band, read band
    # into an array and put array into a list
    if bandnum > 0:
        for files in filelist:

            dataset = gdal.Open(files, GA_ReadOnly)

            band = dataset.GetRasterBand(bandnum)

            bandval = band.ReadAsArray(0,0, band.XSize, band.YSize)

            gribdatalist.append(bandval)

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
    data_dir : string
        Path to data location

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


def windcontours(basemap, uband, vband, latitudes, longitudes, ax=None,
                 figsize=(20,20), wmin=20.0, wmax=85.0, wstep=5.0,
                 colormap='jet', zorder=13, limits=None):
    """
    Contour map of windspeed.

    Parameters
    ----------
    basemap : MapDesign or Basemap instance
        MapDesign instances contain the parameters to initialize a Basemap.
        If a MapDesign is provided, the Basemap may be initialized on a
        given axis using the kwarg `ax`.  If ax is None, then new
        figure, axis, and Basemap instances will be created.
        A previously-generated Basemap may be provided instead of a
        MapDesign.
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
    ax : matplotlib axes instance
        Default None.  The axis on which to draw a new Basemap, if
        `basemap` is not a Basemap instance.  If None, and `basemap` is a
        MapDesign instance, then new figure and axis instances
        will be created.
    figsize : tuple of ints
        Default (20,20).  The size of a new figure instance, if one must be
        created.
    wmin : float
        Default 0.0.  The minimum value for creating contour levels and
        colormapping.  If None, wmin will be the minimum windspeed.
    wmax : float
        Default None.  The maximum value for creating contour levels and
        colormapping.  If None, wmax will be the maximum windspeed.
    wstep : int or float
        The step size between contour levels bounded by wmin, wmax.
    colormap : string
        Default 'jet'.  ['jet'|'blues'|'anomaly'|'heat'|'earth']
        Colormap for windspeed data
    zorder : int
        Default 13.  The zorder of the contours.
    limits : list of floats or ints
        Default None.  Used to trim zdata grid, longitudes, latitudes.
        [bottom latitude, left longitude, top latitude, right longitude]

    Returns
    -------
    fig : matplotlib Figure instance, optional
        Newly created Figure instance.  Only returned if a new figure
        and axis had to be created, i.e. ax=None and `basemap`
        is a MapDesign instance.
    ax : matplotlib Axes instance, optional
        Newly created Axes instance.  ONly returned if a new figure
        and axis had to be created, i.e. ax=None and `basemap`
        is a MapDesign instance.
    windmap : Basemap instance
        Basemap instance with filled contour plot of windspeed data
    w : matplotlib PathCollection instance
        Mappable for use in creating colorbars.  Colorbars may be created
        using make_cbar() or make_cax_cbar().

    """

    # Initialize basemap
    try:
        if ax is None:
            fig, ax, windmap = basemap.make_basemap(figsize)
        else:
            windmap = basemap.make_basemap(figsize, ax=ax)
    except AttributeError:
        windmap = basemap

    # Trim grid to a more localized box
    if limits is not None:
        _, _, uband = gridlimit(limits, longitudes, latitudes, uband)
        longitudes, latitudes, vband = gridlimit(limits, longitudes, latitudes)

    # Calculate windspeed
    windspeed = np.sqrt(uband**2 + vband**2)

    # Inizialize contour levels
    levels = np.arange(wmin, wmax+wstep, wstep)

    # Make grid for contour plot
    lons, lats = np.meshgrid(longitudes, latitudes)

    # Plot data as filled contours
    w = windmap.contourf(lons, lats, windspeed, zorder=contour_zorder,
                         cmap=mm.get_colormap(colormap),
                         levels=levels, latlon=True)

    try:
        return fig, ax, windmap, w
    except:
        return windmap, w


def windbarbs(basemap, uband, vband, latitudes, longitudes, ax=None,
              figsize=(20,20), vectors='arrows', zorder=20):
    """
    Map wind velocity using vectors (arrows or barbs.

    Parameters
    ----------
    basemap : MapDesign or Basemap instance
        MapDesign instances contain the parameters to initialize a Basemap.
        If a MapDesign is provided, the Basemap may be initialized on a
        given axis using the kwarg `ax`.  If ax is None, then new
        figure, axis, and Basemap instances will be created.
        A previously-generated Basemap may be provided instead of a
        MapDesign.
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
    ax : matplotlib axes instance
        Default None.  The axis on which to draw a new Basemap, if
        `basemap` is not a Basemap instance.  If None, and `basemap` is a
        MapDesign instance, then new figure and axis instances
        will be created.
    figsize : tuple of ints
        Default (20,20).  The size of a new figure instance, if one must be
        created.
    vectors : string
        Default 'arrows'.  ['arrows'|'barbs']
        The vector type to indicate wind velocity.  'arrows' comes with a key.
    zorder : int
        Default 20.  The zorder of the vectors.
    limits : list of floats or ints
        Default None.  Used to trim zdata grid, longitudes, latitudes.
        [bottom latitude, left longitude, top latitude, right longitude]

    Returns
    -------
    fig : matplotlib Figure instance, optional
        Newly created Figure instance.  Only returned if a new figure
        and axis had to be created, i.e. ax=None and `basemap`
        is a MapDesign instance.
    ax : matplotlib Axes instance, optional
        Newly created Axes instance.  ONly returned if a new figure
        and axis had to be created, i.e. ax=None and `basemap`
        is a MapDesign instance.
    windmap : Basemap instance
        Basemap instance with vector plot of windspeed data


    """

    # Initialize basemap
    try:
        if ax is None:
            fig, ax, windmap = basemap.make_basemap(figsize)
        else:
            windmap = basemap.make_basemap(figsize, ax=ax)
    except AttributeError:
        windmap = basemap

    # Trim grid to a more localized box
    if limits is not None:
        _, _, uband = gridlimit(limits, longitudes, latitudes, uband)
        longitudes, latitudes, vband = gridlimit(limits, longitudes, latitudes)

    # Change longitudes to -180 to 180
    for i in range(0, longitudes.size):
        if longitudes[i]>180.0:
            longitudes[i] = longitudes[i] - 360.0

    longitude_sort = np.argsort(longitudes)

    longitudes = longitudes[longitude_sort]
    uband = uband[:, longitude_sort]
    vband = vband[:, longitude_sort]


    # # Make grid for wind barbs/ vectors
    # # First, rearrange grid so that it start with a positive longitude
    # firstlon = longitudes[0] + 360.0
    # longitudes = longitudes.tolist()[1:]
    # longitudes.append(firstlon)
    # longitudes = np.asarray(longitudes)

    # firstucol = np.atleast_2d(uband[:,0])
    # uband = np.concatenate((uband[:,1:], firstucol.T), axis=1)

    # firstvcol = np.atleast_2d(vband[:,0])
    # vband = np.concatenate((vband[:,1:], firstvcol.T), axis=1)

    # # Shift grid so longitudes are increasing from -180 to 180
    # # if the projection is not cylindrical
    # if not mapdesign_obj.projection[:-3] is 'cyl':
    #     uband, newlons = shiftgrid(longitudes[longitudes.size/2], uband,
    #                                longitudes, start=False)
    #     vband, newlons = shiftgrid(longitudes[longitudes.size/2], vband,
    #                                longitudes, start=False)
    #     longitudes = newlons

    # Transform vectors to lie correctly on projection
    u, v, lns, lts = windmap.transform_vector(uband[::-1,:], vband[::-1,:],
                                              longitudes, latitudes[::-1],
                                              12, 16, returnxy=True)
    # Plot vectors as arrows or barbs
    if vectors is 'arrows':
        q = windmap.quiver(lns, lts, u, v, zorder=vector_zorder, scale=1000)
        qkey = plt.quiverkey(q, 0.25, 0.05, 60, '60 m/s', labelpos='W',
                             labelcolor='red', fontproperties={'size':20},
                             zorder=zorder, color='red')

    elif vectors is 'barbs':
        q = windmap.barbs(lns, lts, u, v, zorder=zorder)

    try:
        return fig, ax, windmap
    except:
        return windmap


def plot_othervar(basemap, vardata, latitudes, longitudes, ax=None,
                  figsize=(20,20), filetype='ncar', color_min=None,
                  color_max=None, color_levelnum=50, colormap='jet',
                  zorder=13, limits=None):
    """
    For plotting any variable that is not wind.  Interrogates ncar (grib) or
        trmm (hdf) files on any surface and plots data from specified years or
        interval of years in a filled contour plot on a basemap.

    Parameters
    ----------
    basemap : MapDesign or Basemap instance
        MapDesign instances contain the parameters to initialize a Basemap.
        If a MapDesign is provided, the Basemap may be initialized on a
        given axis using the kwarg `ax`.  If ax is None, then new
        figure, axis, and Basemap instances will be created.
        A previously-generated Basemap may be provided instead of a
        MapDesign.
    vardata : (M, N) ndarray of floats
        Data array
    latitudes : (M) ndarray of floats
        Array of latitudes of vardata
    longitudes : (N) ndarray of floats
        Array of longitudes of vardata

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
    filetype : string
        Default 'ncar'.  ['ncar'|'trmm'].  TRMM data has a different
        orientation than NCAR data and must be transposed before plotting.
    color_min : int or float
        Default None.  The minimum value for color mapping.  If None,
        color_min will be the minimum value of the data.
    color_max : int or float
        Default None.  The maximum value for color mapping.  If None,
        color_max will be the maximum value of the data.
    color_levelnum
    colormap : string
        Default 'jet'.  ['jet'|'blues'|'anomaly'|'heat'|'earth']
        Colormap for zdata
    zorder : int
        Default 13.  The zorder of the contours.
    limits : list of floats or ints
        Default None.  Used to trim zdata grid, longitudes, latitudes.
        [bottom latitude, left longitude, top latitude, right longitude]

    Returns
    -------
    fig : matplotlib Figure instance, optional
        Newly created Figure instance.  Only returned if a new figure
        and axis had to be created, i.e. ax=None and `basemap`
        is a MapDesign instance.
    ax : matplotlib Axes instance, optional
        Newly created Axes instance.  ONly returned if a new figure
        and axis had to be created, i.e. ax=None and `basemap`
        is a MapDesign instance.
    varmap : Basemap instance
        Basemap instance will filled contour plot of vardata
    v : matplotlib PathCollection instance
        Mappable for use in creating colorbars.  Colorbars may be created
        using make_cbar() or make_cax_cbar().

    """

    # Initalize basemap
    try:
        if ax is None:
            fig, ax, varmap = basemap.make_basemap(figsize)
        else:
            ax, varmap = basemap.make_basemap(figsize, ax=ax)
    except AttributeError:
        varmap = basemap

    # Trim grid to a more localized box
    if limits is not None:
        longitudes, latitudes, vardata = gridlimit(limits, longitudes,
                                                   latitudes, vardata)

    if color_min is None:
        color_min = varband.min()

    if color_max is None:
        color_max = varband.max()

    levels = np.linspace(color_min, color_max, color_levelnum)

    # Make grid
    lons, lats = np.meshgrid(longitudes, latitudes)

    # TRMM files are oriented weirdly
    if filetype is 'trmm':
        varband = varband.T

    # Make the contour plot
    v = varmap.contourf(lons, lats, varband, cmap=get_colormap(colormap),
                        zorder=zorder, levels=levels, latlon=True)

    try:
        return fig, ax, varmap, v
    except:
        return varmap, v


def zplot(basemap, zdata, latitudes, longitudes, contours, ax=None,
          figsize=(20,20), is_geopotential=True, zmin=0.0, zmax=None,
          zstep=100.0, colormap='jet', colors=None, linewidths=2, zorder=14,
          limits=None):
    """
    Map all or specific z (or any) contours.  Converts geopotential to meters
        if necessary.

    Parameters
    ----------
    basemap : MapDesign or Basemap instance
        MapDesign instances contain the parameters to initialize a Basemap.
        If a MapDesign is provided, the Basemap may be initialized on a
        given axis using the kwarg `ax`.  If ax is None, then new
        figure, axis, and Basemap instances will be created.
        A previously-generated Basemap may be provided instead of a
        MapDesign.
    zdata : (M, N) ndarray of floats
        Data array
    latitudes : (M) ndarray of floats
        Array of latitudes of zdata
    longitudes : (N) ndarray of floats
        Array of longitudes of zdata
    contours : string or list of ints or floats
        'all', or a list of particular contours to plot.

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
    is_geopotential : Boolean
        Default True.  If True, data will be divided by gravitational
        acceleration to get elevation in meters
    zmin : float
        Default 0.0.  The minimum value for creating contour levels and maybe
        colormapping.
    zmax : float
        Default None.  The maximum value for creating contour levels and maybe
        colormapping.  If None, zmax will be the maximum value of zdata.
    zstep : int or float
        The step size between contour levels bounded by zmin, zmax.
    colormap : string
        Default 'jet'.  ['jet'|'blues'|'anomaly'|'heat'|'earth'].  The colormap
        for plotting 'all' contour data.
    colors : list of strings
        Default None.  The colors of particular contours to plot.
        len(colors) == len(contours)
    linewidths : int or float or list of ints or floats
        Default 2.  The linewidth of all contours, or a list of linewidths for
        particularly contours.  Must be in list format if not 'all' contours
        are chosen.  len(linewidths) == len(contours)
    zorder : int
        Default 14.  The zorder of the contour data
    limits : list of floats or ints
        Default None.  Used to trim zdata grid, longitudes, latitudes.
        [bottom latitude, left longitude, top latitude, right longitude]
    """

    # Intialize zmap
    try:
        if ax is None:
            fig, ax, zmap = basemap.make_basemap(figsize)
        else:
            ax, zmap = basemap.make_basemap(figsize, ax=ax)
    except AttributeError:
        zmap = basemap

    # Calculate meters from geopotential
    if is_geopotential:
        zdata = zdata / 9.80665

    # Trim grid to a more localized box
    if limits is not None:
        longitudes, latitudes, zdata = gridlimit(limits, longitudes, latitudes,
                                                 zdata)
    if zmax is None:
        zmax = zdata.max()

    lons, lats = np.meshgrid(longitudes, latitudes)

    levels = np.arange(zmin, zmax, zstep)

    if contours is 'all':
        colormap = mm.get_colormap(colormap)

        z = zmap.contour(lons, lats, zdata, cmap=colormap, zorder=zorder,
                         linewidth=linewidths, latlon=True, levels=levels)

        try:
            return fig, ax, zmap, z
        except:
            return zmap, z

    # PLot specific contours
    else:
        z = zmap.contour(lons, lats, zdata, colors='r', levels=levels,
                         zorder=zorder, latlon=True)

        for lvl in range(0, z.levels.size):
            zl = z.collections[lvl]

            if levels[lvl] in contours:
                ind = contours.index(levels[lvl])
                plt.setp(zl, color=colors[ind], linewidth=linewidths[ind])
            else:
                plt.setp(zl, color=None, alpha=0)

        try:
            return fig, ax, zmap
        except:
            return zmap


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
        lat_inds = np.nonzero(latitudes>limits[0])[0]
        latitudes = latitudes[lat_inds]
        data = data[lat_inds,:]

    if limits[2] is not None:
        lat_inds = np.nonzero(latitudes<limits[2])[0]
        latitudes = latitudes[lat_inds]
        data = data[lat_inds, :]

    if limits[1] > limits[3] and limits[1] > 180:
        lon_tmp = []
        for lon in longitudes:
            if lon> 180:
                lon = 360 - lon
            lon_tmp.append(lon)
        longitudes = np.asarray(lon_tmp)

    if limits[1] is not None:
        lon_inds = np.nonzero(longitudes>limits[1])[0]
        longitudes = longitudes[lon_inds]
        data = data[:, lon_inds]

    if limits[3] is not None:
        lon_inds = np.nonzero(longitudes<limits[3])[0]
        longitudes = longitudes[lon_inds]
        data = data[:, lon_inds]

    return longitudes, latitudes, data