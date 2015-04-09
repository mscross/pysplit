from __future__ import division
import numpy as np
import grib_reader as gr


def oni_file_reader(oni_path):
    """
    Reads the Oceanic Nino Index from a text file into a NumPy array.

    Parameters
    ----------
    oni_path : string
        File name and location of ONI data file in plaintext.

    Returns
    -------
    year_list : list of ints
        Years with complete data in ONI file
    oni_array : 2D Numpy ndarray of floats
        ONI data arranged by year (row) and three month period (column)
    trimonth_list : list of strings
        Column header for oni_array.  List of three character strings ('DJF')
        representing three month periods.  Length = 12

    Notes
    -----
    Oceanic Nino Index data may be downloaded `here`_:

    .. _here: http://www.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ensoyears.shtml

    Data must be saved in a plaintext file before use.

    """

    year_list = []
    oni_array = np.empty((0, 12))

    # Open file
    with open(oni_path, 'r') as oni_file:

        # Get header
        headerline = oni_file.readline()
        trimonth_list = headerline.split()[1:]

        # Read through rest of file until at end or incomplete year encountered
        while True:
            oniline = oni_file.readline()

            # Break if at end of file
            if oniline == '':
                break

            data = oniline.split()

            # Break if an incomplete year is encountered
            if len(data) < 13:
                break

            # Put year into list
            year_list.append(data[0])

            # Put data into array
            onidata = np.asarray(data[1:]).astype(np.float64)
            onidata = np.atleast_2d(onidata)
            oni_array = np.concatenate((oni_array, onidata), axis=0)

    # Change list of strings to list of ints
    year_list = [int(i) for i in year_list]

    return year_list, oni_array, trimonth_list


def oni_file_interrogator(monthstring, phase, strength, oni_array,
                          year_list, trimonth_list):
    """
    Find years a given three-month period is of the ENSO phase & strength of
    interest.

    Parameters
    ----------
    monthstring : string
        Three letter string representing three contiguous months
    phase : string
        The phase of interest.  ['el nino'|'la nina'|'none']
    strength : string
        Ignored if `phase` is 'none'.
        ['strong'|'moderate'|'weak'|'all'|'high'|'low'].  'high' and 'low'
        encompass 'strong' and 'moderate', 'moderate' and 'weak'
    oni_array : 2D Numpy ndarray of floats
        ONI data [year, 3-month period]
    year_list : list of ints
        Years with complete data in ONI file
    trimonth_list : list of strings
        Column header for oni_array.  List of three character strings ('DJF')
        representing three month periods.  Length = 12

    Returns
    -------
    years_ofinterest : Numpy ndarray of ints
        Years where interesting period is the phase and strength(s) of interest

    """

    # Dictionary of upper and lower limits for each phase and strength
    # or range of strength
    enso_dict = {'el nino': {'strong': (1.5, 3.5),
                             'moderate': (1.0, 1.4),
                             'weak': (0.5, 0.9),
                             'all': (0.5, 3.5),
                             'high': (1.0, 3.5),
                             'low': (0.5, 1.4)},
                 'la nina': {'strong': (-3.5, -1.5),
                             'moderate': (-1.4, -1.0),
                             'weak': (-0.9, -0.5),
                             'all': (-3.5, -0.5),
                             'high': (-3.5, -1.0),
                             'low': (-1.4, -0.5)},
                 'none' : (-0.4, 0.4)}

    # Get the column number of the data corresponding to the 3 month period
    trimonth_index = trimonth_list.index(monthstring)

    # Get the upper and lower limits of the phase and strength of interest
    try:
        limits = enso_dict[phase][strength]
    except:
        limits = enso_dict[phase]

    # Get the row indices of the appropriate column where the ONI value is
    # within the upper and lower limits determined above
    yargs = np.nonzero((oni_array[:, trimonth_index] >= limits[0]) &
                       (oni_array[:, trimonth_index] <= limits[1]))[0]

    # Turn the given list of years with data into an array of ints
    year_array = np.asarray(year_list).round().astype(int)

    # Get an array of ints where each item is only one of the years
    # that meet the phase and strength/ strength range criteria
    years_ofinterest = year_array[yargs]

    return years_ofinterest


def find_years(oni_path, monthstring, phase, strength,
               data_lowerlimit=1980, data_upperlimit=2012):
    """
    Takes the phase & phase strength and the three month period of interest,
    gets the array of Oceanic Nino Index values via ``oni_file_reader()``,
    gets the years of interest via ``oni_file_interrogator()``, limits those
    years to ranges covered by the chosen dataset, and returns the list of
    single-item lists of the 3 month numbers of interest, plus a list of lists
    of years of interest corresponding to each month (monthstring = 'DJF', for
    example, has December of the previous year)

    Parameters
    ----------
    oni_path : string
        file name and location of ONI data text file
    monthstring : string
        Three letter string representing three contiguous months
    phase : string
        The phase of interest.  ['el nino'|'la nina'|'none']
    strength : string
        Ignored if `phase` is 'none'.
        ['strong'|'moderate'|'weak'|'all'|'high'|'low'].  'high' and 'low'
        encompass 'strong' and 'moderate', 'moderate' and 'weak'

    Keyword Arguments
    -----------------
    data_lowerlimit : int
        Default 1979.  The oldest year of ENSO information requested/available
    data_upperlimit : int
        Default 2012.  The most recent year of ENSO data requested/available

    Returns
    -------
    month_list : list of lists of 1 int each
        ``gr.file_dates()`` requires months in a list format.
        Nested listing permits each month to be passed separately, which is
        necessary when the three month period of interest is one that spans
        the turn of the year.
    years_list : list of lists of ints
        Each sub-list contains the years where the phase/strength/data limit
        conditions were met for one of the months in the three month period of
        interest

    """

    # Get array of Oceanic Nino Index data, plus the list of years (rows)
    # and list of three-month-groups (columns, 12)
    year_list, oni_array, trimonth_list = oni_file_reader(oni_path)

    # Initialize dictionary
    monthlist_dict = {'DJF': [[12], [1],  [2]],  'JFM': [[1],  [2],  [3]],
                      'FMA': [[2],  [3],  [4]],  'MAM': [[3],  [4],  [5]],
                      'AMJ': [[4],  [5],  [6]],  'MJJ': [[5],  [6],  [7]],
                      'JJA': [[6],  [7],  [8]],  'JAS': [[7],  [8],  [9]],
                      'ASO': [[8],  [9],  [10]], 'SON': [[9],  [10], [11]],
                      'OND': [[10], [11], [12]], 'NDJ': [[11], [12], [1]]}

    singlemonth_dict = {'1' : 'DJF', '2' : 'JFM', '3' : 'FMA', '4' : 'MAM',
                        '5' : 'AMJ', '6' : 'MJJ', '7' : 'JJA', '8' : 'JAS',
                        '9' : 'ASO', '10': 'SON', '11': 'OND', '12': 'NDJ'}

    # Convert user input to appropriate case
    phase = phase.lower()
    strength = strength.lower()

    # Get list of lists of months
    if len(monthstring) > 1:
        threemonths = True
        month_list = monthlist_dict[monthstring]
    else:
        threemonths = False
        month_list = [[int(monthstring)]]
        monthstring = singlemonth_dict[monthstring]

    # Get the list of years where the phase/strength criteria is met
    years_ofinterest = oni_file_interrogator(monthstring, phase, strength,
                                             oni_array, year_list,
                                             trimonth_list)

    years_ofinterest_adj = []

    # Limit data to a particular time period
    years_ofinterest = years_ofinterest[years_ofinterest >= data_lowerlimit]
    years_ofinterest = years_ofinterest[years_ofinterest <= data_upperlimit]

    # Turn array of ints into a list
    years_ofinterest = years_ofinterest.tolist()

    if threemonths:
        # Two of three month periods include a month from previous or next year
        # Construct a list of the lists of years
        if monthstring is 'DJF':
            years_ofinterest_adj = [i - 1 for i in years_ofinterest]
            years_list = [years_ofinterest_adj, years_ofinterest,
                          years_ofinterest]

        elif monthstring is 'NDJ':
            years_ofinterest_adj = [i + 1 for i in years_ofinterest]
            years_list = [years_ofinterest, years_ofinterest,
                          years_ofinterest_adj]

        else:
            years_list = [years_ofinterest, years_ofinterest, years_ofinterest]
    else:
        years_list = [years_ofinterest]

    # print years_list

    return month_list, years_list


def enso_winddata(month_list, years_list, level, data_dir, startstring):
    """
    Prepare for wind plot by gathering data representing the ENSO phase
        and strength of interest

    Can find anomaly by performing file_years(), enso_winddata() twice with
        different phase, strength

    Parameters
    ----------
    month_list : list of 1 or 3 lists of 1 int
        Each month from a one or three-month period in a separate list
    years_list : list of 1 or 3 lists of ints
        The years that had the interesting phase/strength during each of the
        months in list.
    level : string
        Desired surface level of data
    data_dir : string
        Full or relative path to the data directory.
    startstring : string
        Initial part of filenames

    Returns
    -------
    ubanddata : (M, N) ndarray of floats
        The U component of the wind data
    vbanddata : (M, N) ndarray of floats
        The V component of the wind data
    lats : (N) ndarray of floats
        The latitudes of the data grid
    lons : (M) ndarray of floats
        The longitudes of the data grid

    """

    ubandlist = []
    vbandlist = []

    for monlist, yrlist in zip(month_list, years_list):

        dates = gr.file_dates(yrlist, False, monlist, False)

        uband, lats, lons = gr.get_gribdata(level, 'U', startstring, data_dir,
                                            dates)
        vband, _, _ = gr.get_gribdata(level, 'V', startstring, data_dir, dates)

        ubandlist.extend(uband)
        vbandlist.extend(vband)

    ubanddata = gr.averager(ubandlist)
    vbanddata = gr.averager(vbandlist)

    return ubanddata, vbanddata, lats, lons


def enso_vardata(month_list, years_list, var, level, filetype,
                 data_dir, startstring):
    """
    Prepare for plot by gathering data representing the ENSO phase
        and strength of interest

    Can find anomaly by performing file_years(), enso_vardata() twice with
        different phase, strength

    Parameters
    ----------
    month_list : list of 1 or 3 lists of 1 int
        Each month from a one or three-month period in a separate list
    years_list : list of 1 or 3 lists of ints
        The years that had the interesting phase/strength during each of the
        months in list.
    var : string
        The variable to inspect
    level : string
        Desired surface level of data
    data_dir : string
        Full or relative path to the data directory.
    startstring : string
        Initial part of filenames

    Returns
    -------
    varbanddata : (M, N) ndarray of floats
        The data
    lats : (N) ndarray of floats
        The latitudes of the data grid
    lons : (M) ndarray of floats
        The longitudes of the data grid

    """

    vardatalist = []

    for monlist, yrlist in zip(month_list, years_list):

        dates = gr.file_dates(yrlist, False, monlist, False)

        if filetype is 'HDF':
            varlist, lons, lats = gr.get_hdfdata(dates, startstring, data_dir,
                                                 var=var)
        else:
            varlist, lats, lons = gr.get_gribdata(level, var, startstring,
                                                  data_dir, dates)

        vardatalist.extend(varlist)

    varbanddata = gr.averager(vardatalist)

    return varbanddata, lats, lons
