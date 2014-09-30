import numpy as np
import math
import grib_reader as gr


def oni_file_reader(working_dir):
    """
    Reads the Oceanic Nino Index from a text file into a Numpy array.

    Parameters
    ----------
    working_dir : string
        file name and location of ONI data text file

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
    Oceanic Nino Index data may be downloaded here:
        http://www.cpc.ncep.noaa.gov/products/
                analysis_monitoring/ensostuff/ensoyears.shtml
    Data must be saved in a text file before use.

    """

    #Initialize
    year_list = []
    oni_array = np.empty((0, 12))

    # Open file
    oni_file = open(working_dir, 'r')

    # Get header
    headerline = oni_file.readline()
    trimonth_list = headerline.split()[1:]

    # Read through rest of file until the file is done
    # or an incomplete year is encountered
    while True:
        oniline = oni_file.readline()

        # Break if at end of file
        if oniline =='':
            break

        data = oniline.split()

        # Break if an incomplete year is encountered
        if len(data)< 13:
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
    Takes an array of data and indexing information and returns all
        the years that the given three month period is of the ENSO phase
        and strength of interest.

    Parameters
    ----------
    monthstring : string
        Three letter string representing three contiguous months
    phase : string
        'el nino', 'la nina' or 'none'.  The phase of interest
    strength : string
        'strong', 'moderate', 'weak', 'all', 'high', 'low'.  Ignored if phase
        is 'none'. 'all', 'high', and 'low' indicate all, strong & moderate,
        and moderate & low strengths
    oni_array : 2D Numpy ndarray of floats
        ONI data arranged by year (row) and three month period (column)
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
    if phase != 'none':
        limits = enso_dict[phase][strength]
    else:
        limits = enso_dict[phase]

    # Get the row indices of the appropriate column where the ONI value is
    # within the upper and lower limits determined above
    yargs = np.nonzero((oni_array[:, trimonth_index] >= limits[0]) &
                       (oni_array[:, trimonth_index] <= limits[1]))[0]

    # Turn the given list of years with data into an array of ints
    year_array = np.asarray(year_list).astype(int)

    # Get an array of ints where each item is only one of the years
    # that meet the phase and strength/ strength range criteria
    years_ofinterest = year_array[yargs]

    return years_ofinterest


def enso_plotprep(working_dir, monthstring, phase, strength,
                  data_lowerlimit=1979, data_upperlimit=2012):
    """
    Takes the phase & phase strength and the three month period of interest,
        gets the array of Oceanic Nino Index values via oni_file_reader(),
        gets the years of interest via oni_file_interrogator(),
        limits those years to years that are covered by the chosen
        dataset, and returns the list of single-item lists of the 3 month
        numbers of interest, plus a list of lists of years of interest
        corresponding to each month (monthstring = 'DJF', for example,
        has December of the previous year)

    Parameters
    ----------
    working_dir : string
        file name and location of ONI data text file
    monthstring : string
        Three letter string representing three contiguous months or
        1 character string representing a single month
    phase : string
        'el nino', 'la nina' or 'none'.  The phase of interest
    strength : string
        'strong', 'moderate', 'weak', 'all', 'high', 'low'.
        Ignored if phase is 'none'.
        'all', 'high', and 'low' indicate all, strong & moderate,
        and moderate & low strengths
    data_lowerlimit : int
        Default 1979.  The oldest year of ENSO information required,
        either because of data limits or user is interested
        in a particular time period
    data_upperlimit : int
        Default 2012.  The most recent year of ENSO information required,
        either because of data limits or user is interested
        in a particular time period

    Returns
    -------
    monthlist : list of lists of 1 int each
        gr.grib_lister() and gr.hdf_lister() require months in a list format.
        Nested listing permits each month to be passed separately, which is necessary
        when the three month period of interest is one that spans the turn of the year.
    years_list : list of lists of ints
        Each sub-list contains the years where the phase/strength/data limit conditions
        were met for one of the months in the three month period of interest

    """

    # Get array of Oceanic Nino Index data, plus the list of years (rows)
    # and list of three-month-groups (columns, 12)
    year_list, oni_array, trimonth_list = oni_file_reader(working_dir)

    # Initialize dictionary
    monthlist_dict = {'DJF': [[12],[1], [2]], 'JFM': [[1], [2], [3]],
                      'FMA': [[2], [3], [4]], 'MAM': [[3], [4], [5]],
                      'AMJ': [[4], [5], [6]], 'MJJ': [[5], [6], [7]],
                      'JJA': [[6], [7], [8]], 'JAS': [[7], [8], [9]],
                      'ASO': [[8], [9], [10]], 'SON': [[9], [10],[11]],
                      'OND': [[10],[11],[12]], 'NDJ': [[11],[12],[1]]}

    singlemonth_dict = {'1' : 'DJF', '2' : 'JFM', '3' : 'FMA', '4' : 'MAM',
                        '5' : 'AMJ', '6' : 'MJJ', '7' : 'JJA', '8' : 'JAS',
                        '9' : 'ASO', '10': 'SON', '11': 'OND', '12': 'NDJ'}

    # Convert user input to appropriate case
    phase = phase.lower()
    strength = strength.lower()

    # Get list of lists of months
    if len(monthstring) > 1:
        threemonths = True
        monthlist = monthlist_dict[monthstring]
    else:
        threemonths = False
        monthlist = [[int(monthstring)]]
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
        # Two of the three month periods include a month from the previous or next year
        # Construct a list of the lists of years
        if monthstring is 'DJF':
            years_ofinterest_adj = [i - 1 for i in years_ofinterest]
            years_list = [years_ofinterest_adj, years_ofinterest, years_ofinterest]

        elif monthstring is 'NDJ':
            years_ofinterest_adj = [i + 1 for i in years_ofinterest]
            years_list = [years_ofinterest, years_ofinterest, years_ofinterest_adj]

        else:
            years_list = [years_ofinterest, years_ofinterest, years_ofinterest]
    else:
        years_list = [years_ofinterest]

    print years_list

    return monthlist, years_list


def enso_winddata(monthstring, phase, strength, level,
                  enso_dir, data_dir,
                  data_lowerlimit=1979, data_upperlimit=2012):
    """
    Prepare for wind plot by gathering data representing the ENSO phase
        and strength of interest

    Parameters
    ----------
    monthstring : string
        Three letter string representing three contiguous months or
        1 character string representing a single month
    phase : string
        'el nino', 'la nina' or 'none'.  The phase of interest
    strength : string
        'strong', 'moderate', 'weak', 'all', 'high', 'low'.
        Ignored if phase is 'none'.
        'all', 'high', and 'low' indicate all, strong & moderate,
        and moderate & low strengths
    level
    enso_dir
    data_dir

    Keyword Arguments
    -----------------
    pl
    measured
    data_lowerlimit : int
        Default 1979.  The oldest year of ENSO information required,
        either because of data limits or user is interested
        in a particular time period
    data_upperlimit : int
        Default 2012.  The most recent year of ENSO information required,
        either because of data limits or user is interested
        in a particular time period

    """

    monthlist, years_list = enso_plotprep(enso_dir, monthstring,
                                          phase, strength, data_lowerlimit,
                                          data_upperlimit)

    ubandlist = []
    vbandlist = []

    for mon, yr, in zip(monthlist, years_list):

        uband, lats, lons = gr.grib_lister(level, 'U', yr, mon, pl,
                                             wind, data_dir, True)
        vband, lats, lons = gr.grib_lister(level, 'V', yr, mon, pl,
                                             wind, data_dir, True)

        ubandlist.extend(uband)
        vbandlist.extend(vband)


    ubanddata = np.mean(np.dstack(ubandlist), axis=2)
    vbanddata = np.mean(np.dstack(vbandlist), axis=2)

    return ubanddata, vbanddata, lats, lons



def enso_vardata(monthstring, phase, strength, var, level, filetype,
                 enso_dir, data_dir,
                 data_lowerlimit=1980, data_upperlimit=2012):
    """
    Prepare for plot by gathering data representing the ENSO phase
        and strength of interest

    Parameters
    ----------
    monthstring : string
        Three letter string representing three contiguous months or
        1 character string representing a single month
    phase : string
        'el nino', 'la nina' or 'none'.  The phase of interest
    strength : string
        'strong', 'moderate', 'weak', 'all', 'high', 'low'.
        Ignored if phase is 'none'.
        'all', 'high', and 'low' indicate all, strong & moderate,
        and moderate & low strengths
    var
    level
    filetype
    enso_dir
    data_dir

    Keyword Arguments
    -----------------
    pl
    measured
    data_lowerlimit : int
        Default 1979.  The oldest year of ENSO information required,
        either because of data limits or user is interested
        in a particular time period
    data_upperlimit : int
        Default 2012.  The most recent year of ENSO information required,
        either because of data limits or user is interested
        in a particular time period

    """

    monthlist, years_list = enso_plotprep(enso_dir, monthstring, phase,
                                          strength, data_lowerlimit,
                                          data_upperlimit)

    varlist = []

    for mon, yr in zip(monthlist, years_list):

        if filetype is 'ncar':
            varband, lats, lons = gr.grib_lister(level, var, yr,mon, pl,
                                                 False, data_dir, measured)

        else:
            varband, lats, lons = gr.hdf_lister(yr, mon, data_dir)

        varlist.extend(varband)

    varbanddata = np.mean(np.dstack(varband), axis=2)


    return varbanddata, lats, lons


def enso_windanomaly(monthstring, phase1, strength1, level, enso_dir,
                     data_dir, phase2='none', strength2='none', pl=True,
                     data_lowerlimit=1980, data_upperlimit=2012):

    """
    Prepare for wind plot by gathering data representing the ENSO phase
        and strength of interest

    Parameters
    ----------
    monthstring : string
        Three letter string representing three contiguous months or
        1 character string representing a single month
    phase1 : string
        'el nino', 'la nina' or 'none'.  The phase of interest
    strength1 : string
        'strong', 'moderate', 'weak', 'all', 'high', 'low'.
        Ignored if phase is 'none'.
        'all', 'high', and 'low' indicate all, strong & moderate,
        and moderate & low strengths
    level
    enso_dir
    data_dir

    Keyword Arguments
    -----------------
    phase2
    strength2
    pl
    measured
    data_lowerlimit : int
        Default 1979.  The oldest year of ENSO information required,
        either because of data limits or user is interested
        in a particular time period
    data_upperlimit : int
        Default 2012.  The most recent year of ENSO information required,
        either because of data limits or user is interested
        in a particular time period

    """

    u1, v1, lats, lons = enso_winddata(monthstring, phase1, strength1, level,
                                       enso_dir, data_dir, pl=pl,
                                       measured=measured,
                                       data_lowerlimit=data_lowerlimit,
                                       data_upperlimit=data_upperlimit)

    u2, v2, lats, lons = enso_winddata(monthstring, phase2, strength2, level,
                                       enso_dir, data_dir, pl=pl,
                                       measured=measured,
                                       data_lowerlimit=data_lowerlimit,
                                       data_upperlimit=data_upperlimit)

    u_anomaly = u1 - u2
    v_anomaly = v1 - v2

    return u_anomaly, v_anomaly, lats, lons


def enso_varanomaly(monthstring, phase1, strength1, var, level, filetype,
                    enso_dir, data_dir, phase2='none', strength2='none',
                    pl=True, measured=True,
                    data_lowerlimit=1980, data_upperlimit=2012):

    """
    Prepare for plot by gathering data representing the ENSO phase
        and strength of interest

    Parameters
    ----------
    monthstring : string
        Three letter string representing three contiguous months or
        1 character string representing a single month
    phase : string
        'el nino', 'la nina' or 'none'.  The phase of interest
    strength : string
        'strong', 'moderate', 'weak', 'all', 'high', 'low'.
        Ignored if phase is 'none'.
        'all', 'high', and 'low' indicate all, strong & moderate,
        and moderate & low strengths
    var
    level
    filetype
    enso_dir
    data_dir

    Keyword Arguments
    -----------------
    phase2
    strength2
    pl
    measured
    data_lowerlimit : int
        Default 1979.  The oldest year of ENSO information required,
        either because of data limits or user is interested
        in a particular time period
    data_upperlimit : int
        Default 2012.  The most recent year of ENSO information required,
        either because of data limits or user is interested
        in a particular time period

    """

    var1, lats, lons = enso_vardata(monthstring, phase1, strength1, var, level,
                                    filetype, enso_dir, data_dir,
                                    pl=pl, measured=measured,
                                    data_lowerlimit=data_lowerlimt,
                                    data_upperlimit=data_upperlimit)

    var2, lats, lons = enso_vardata(monthstring, phase2, strength2, var, level,
                                    filetype, enso_dir, data_dir,
                                    pl=pl, measured=measured,
                                    data_lowerlimit=data_lowerlimt,
                                    data_upperlimit=data_upperlimit)

    var_anomaly = var1 - var2

    return var_anomaly, lats, lons