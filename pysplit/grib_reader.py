import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
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
    filelist : list of files
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
        LIst of strings representing dates in the format YYYYMM

    """

    dates = []

    if years_isrange:
        years = range(years[0], years[-1]+1)
    if months_isrange:
        months = range(months[0], months[-1]+1)

    for y, m in itertools.product(years, months):
        dates.append('{0:04}{1:02}'.format(y, m))

    return dates



def grib_lister(level, var, startstring, data_dir, filedates=None,
                filenames=None):
    """
    Gets a list of GRIB files, figures out which band to look at, then opens
        files and reads in data to an array to be averaged.

    Parameters
    ----------
    level : string
        Desired surface level of data
    var : string
        Band variable
    years : list of ints
        The years of data desired
    months : list of ints
        The month or months of data to include
    data_dir : string
        Path to data location

    Returns
    -------
    gribdata : list of (M, N) ndarray of floats
        Arrays of data from all GRIB files in indicated time period
    latitudes : ndarray of floats
        (M) array of latitudes of gribdata from 89 to -89 (256)
    longitudes : ndarray of floats
        (N) array of longitudes of gribdata from -0.3 to 359.7 (512)

    """

    # Get a list of files that meet the timeframe and data criteria
    if filedates is not None:
        filelist = file_lister(filedates, startstring, data_dir)
    elif filenames is not None:
        filelist = filenames
    else:
        raise ValueError('No data files indicated!')

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


def hdf_lister(filedates, startstring, data_dir):
    """
    Gets a list of HDF files, opens them, reads data into array to be averaged.

    Parameters
    ----------
    years : list of ints
        The years of data to include
    months : list of ints
        The months within the years of data to include
    data_dir : string
        Path to data location

    Returns
    -------
    hdfdata : list of (M, N) ndarray of floats
        arrays of data from all trmm (HDF) files given
    longitudes : (M) ndarray of floats
        array of longitudes of hdfdata from 180 to -180 (1440)
    latitudes : (N) ndarray of floats
        array of latitudes of hdfdata from -50 to 50 (400)

    """

    filelist = file_lister(filedates, startstring, data_dir)
    hdfdatalist = []

    for files in filelist:

        hdf_file = pyhdf.SD.SD(files)

        precip = hdf_file.select('precipitation')

        precipdata = precip.get()

        hdfdatalist.append(precipdata)

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



def windplot(basemap, uband, vband, latitudes, longitudes, figsize=(20,20),
             vectors_only=False, vectors='arrows', scale=(20.0, 85.0, 5.0),
             colormap='jet', contour_zorder=13, vector_zorder=20,
             cbar_size=(20,1.0), reverse_cbar=False, cbar_position='right'):
    """
    Map winds at a particular pressure level.

    Only interrogates GRIB files. Choice of windbarbs or arrows with scale.
        Calls grib_lister to average U and V components over the specified
        period, calls ps.make_basemap to make basemap

    Parameters
    ----------
    basemap : MapDesign object or Basemap object
        Object containing the parameters to intialize Basemap or a
        previously-made Basemap, probably with data already plotted on it
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
    figsize : tuple of ints
        Default (20,20).  Dimensions of figure
    vectors_only : Boolean
        Default False.  If True, wind contour plotting is skipped
    vectors : string
        Default 'arrows'.  Determines whether vectors ('arrows') with a scale or
        wind barbs ('barbs') are plotted
    scale : tuple of floats
        Default (20.0, 85.0, 5.0).  (min, max, step) for contouring
    colormap : string
        Default 'jet'.  ['jet'|'blues'|'anomaly'|'heat'|'earth']
        Passed to a dictionary which retrieves the indicated colormap
    contour_zorder : int
        Default 13.  The zorder of the filled contours.
    vector_zorder : int
        Default 20.  The zorder of the vectors.
    cbar_size : tuple of int, float
        Default (20, 1.0).  Colobar (aspect, shrink).  The H/W ratio of the
        colorbar, the fractional size of the colorbar.
    reverse_cbar : Boolean
        Default False.  If True, colorbar will be reversed (values will
        still map to same colors)
    cbar_position : string
        Default 'right'.  ['bottom'|'right'].  The location of the
        colorbar relative to the map.

    Returns
    -------
    windmap : plot
        Filled contour map of windspeed plus vector plot of wind velocity
        on a Basemap instance

    """

    # Initialize basemap
    try:
        windmap = basemap.make_basemap(figsize)
    except AttributeError:
        windmap = basemap

    # Calculate windspeed
    windspeed = np.sqrt(uband**2 + vband**2)

    # Inizialize contour levels
    levels = np.arange(scale[0], scale[1]+scale[2], scale[2])

    # Make grid for contour plot
    lons, lats = np.meshgrid(longitudes, latitudes)

    if not vectors_only:

        # Plot data as filled contours
        w = windmap.contourf(lons, lats, windspeed, zorder=contour_zorder,
                             cmap=mm.get_colormap(colormap),
                             levels=levels, latlon=True)

        # Set up colorbar
        cbar = windmap.colorbar(w, location=cbar_position, pad='3%',
                                aspect=cbar_size[0], shrink=cbar_size[1])
        cbar.ax.tick_params(labelsize=16)

        mm.edit_cbar(cbar, cbar_position, reverse_cbar, 'Wind Speed (m/s/)')

    # Make grid and plot wind vectors
    if vectors is not None:

        # Make grid for wind barbs/ vectors
        # First, rearrange grid so that it start with a positive longitude
        firstlon = longitudes[0] + 360.0
        longitudes = longitudes.tolist()[1:]
        longitudes.append(firstlon)
        longitudes = np.asarray(longitudes)

        firstucol = np.atleast_2d(uband[:,0])
        uband = np.concatenate((uband[:,1:], firstucol.T), axis=1)

        firstvcol = np.atleast_2d(vband[:,0])
        vband = np.concatenate((vband[:,1:], firstvcol.T), axis=1)

        # Shift grid so longitudes are increasing from -180 to 180
        # if the projection is not cylindrical
        if not mapdesign_obj.projection[:-3] is 'cyl':
            uband, newlons = shiftgrid(longitudes[longitudes.size/2], uband,
                                       longitudes, start=False)
            vband, newlons = shiftgrid(longitudes[longitudes.size/2], vband,
                                       longitudes, start=False)
            longitudes = newlons

        # Transform vectors to lie correctly on projection
        u, v, lns, lts = windmap.transform_vector(uband[::-1,:], vband[::-1,:],
                                                  longitudes, latitudes[::-1],
                                                  12, 16, returnxy=True)
        # Plot vectors as arrows or barbs
        if vectors is 'arrows':
            q = windmap.quiver(lns, lts, u, v, zorder=vector_zorder, scale=1000)
            qkey = plt.quiverkey(q, 0.25, 0.05, 60, '60 m/s', labelpos='W',
                                 labelcolor='red', fontproperties={'size':20},
                                 zorder=vector_zorder, color='red')

        elif vectors is 'barbs':
            q = windmap.barbs(lns, lts, u, v, zorder=vector_zorder)

    return windmap


def plot_othervar(basemap, vardata, latitudes, longitudes,
                  filetype='ncar', figsize=(20,20),
                  scale=(0.0, None, 50), colormap='jet',
                  zorder=13, cbar_label=None, cbar_size=(20,1.0),
                  reverse_cbar=False, cbar_position='right'):
    """
    For plotting any variable that is not wind.  Interrogates ncar (grib) or
        trmm (hdf) files on any surface and plots data from specified years or
        interval of years in a filled contour plot on a basemap.

    Parameters
    ----------
    basemap : MapDesign object or Basemap object
        Object containing the parameters to intialize Basemap or a
        previously-made Basemap, probably with data already plotted on it
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
    figsize : tuple of ints
        Default (20,20).  Dimensions of figure
    scale : tuple of floats
        Default (0.0, None, 50).  (minimum, maximum, # of scale divisions)
    colormap : string
        Default 'jet'.  ['jet'|'blues'|'anomaly'|'heat'|'earth']
        Passed to a dictionary which retrieves the indicated colormap
    zorder : int
        Default 13.  The zorder of the contours.
    cbar_label : string
        Default None.  Label for colorbar
    cbar_size : tuple of int, float
        Default (20, 1.0).  Colobar (aspect, shrink).  The H/W ratio of the
        colorbar, the fractional size of the colorbar.
    reverse_cbar : Boolean
        Default False.  If True, colorbar will be reversed (values will
        still map to same colors)
    cbar_position : string
        Default 'bottom'.  ['bottom'|'right'].  The location of the
        colorbar relative to the map.

    Returns
    -------
    varmap : plot
        Filled contour plot of user-specified variable on a basemap instance

    """

    # Initalize basemap
    try:
        varmap = basemap.make_basemap(figsize)
    except AttributeError:
        varmap = basemap

    # Set the contour levels
    if scale[1] is None:
        scale[1] = varband.max()

    levels = np.linspace(scale[0], scale[1], scale[2])

    # Make grid
    lons, lats = np.meshgrid(longitudes, latitudes)

    # TRMM files are oriented weirdly
    if filetype is 'trmm':
        varband = varband.T

    # Make the contour plot
    v = varmap.contourf(lons, lats, varband,
                        cmap=ps.get_colormap(colormap),
                        zorder=zorder, levels=levels, latlon=True)

    # Set up the colorbar and ticks
    cbar = varmap.colorbar(v, location=cbar_position, pad='5%',
                           aspect=cbar_size[0], shrink=cbar_size[1])

    cbar.locator = tk.MaxNLocator(ticks, integer=False)
    cbar.ax.tick_params(labelsize=16)
    cbar.update_ticks()

    # Edit the colorbar direction and label orientation
    mm.edit_cbar(cbar, cbar_position, reverse_cbar, cbar_label)

    return varmap