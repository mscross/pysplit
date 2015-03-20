import os
from subprocess import call
import itertools


def generate_trajectories(basename, hysplit_working, output_dir, meteo_path,
                          years, years_isrange, months, months_isrange, hours,
                          altitudes, coordinates, run, isbackward,
                          meteo_type='gdas1', get_forward=True,
                          get_clippedtraj=True):
    """
    Run sequence of HYSPLIT simulations over a given time and set of altitudes.

    This does not set along-trajectory meteorological output- edit `SETUP.CFG`
        in the HYSPLIT working directory or in the HYSPLIT4 GUI to reflect
        desired output variables.  It is also recommended to change TRATIO
        in the SETUP file to 0.25 to limit integration error.

    To estimate integration error for back trajectory calculations, forward
        trajectories must be generated from the start of back trajectories
        during calculation of back trajectories.

    One trajectory per simulation file is output.  Output may be multiline,
        depending on number of variables selected in SETUP (maximum 7 of 9
        may be selected for single-line output).  Multiline output not
        supported in some applications (clustering).  To cluster multiline
        output, set get_clippedtraj=True.  Clipped trajectories have all data
        removed except for their path information, so they are single-line
        output and supported for clustering operations.

    Parameters
    ----------
    basename : string
        Basename for all output files.
    hysplit_working : string
        Full or relative path to the HYSPLIT working directory.
    output_dir : string
        Full or relative path to the desired output directory.
    meteo_path : string
        Full or relative path to the location of the meteorology files
    years : list of ints
        The year(s) to run simulations.  Can be a list of specific years
        (set `years_isrange` to False) or the start and end of an inclusive
        range of years (set `years_isrange` to True)
    years_isrange : Boolean
        Indicates if `years` is a list to be used as is or is a range to
        be generated.
    months : list of ints
        The month(s) to run simulations.  Can be a list of specific months
        (set `months_isrange` to False) or the start and end of an inclusive
        range of months (set `months_isrange` to True)
    months_isrange : Boolean
        Indicates if `months` is a list to be used as is or a range to
        be generated.
    hours : list of ints
        Parcel launching times in UTC.
    altitudes : list of ints
        The altitudes in meters above ground level to launch parcels from
    coordinates : tuple of floats
        The parcel starting location in decimal degrees.
        Format is (latitude, longitude)
    run : int
        Length (hours) of simulation.  If necessary, sign is corrected to
        match isbackward argument.
    isbackward : Boolean
        Indicates back trajectory calculation

    Keyword Arguments
    -----------------
    meteo_type : string
        Default 'gdas1'.  The type of meteorology to use.
    get_forward : Boolean
        Default True.  [True|False].  If True and isbackward is also True,
        then a forward trajectory will be calculated from earliest point of
        each back trajectory and stored in a subfolder in output_dir
    get_clippedtraj : Boolean
        Default True.  [True|False].  Outputs trajectory files with
        single-line timesteps containing only path information.  Provides
        clustering support to multiline files.

    Returns
    -------
    Nothing, HYSPLIT is called and results are written to `output_dir`.

    """

    # Hardcode "CONTROL" as this is required
    filename = "CONTROL"
    orig_dir = os.getcwd()

    # Simulation direction check
    if isbackward:
        if run > 0:
            run = run * -1
    else:
        if run < 0:
            run = run * -1
        get_forward = False

    if years_isrange:
        years = range(years[0], years[-1] + 1)

    if months_isrange:
        months = range(months[0], months[-1] + 1)

    season_month_days = {12: ['winter', 'dec', 31],
                         1 : ['winter', 'jan', 31],
                         2 : ['winter', 'feb', 28],
                         3 : ['spring', 'mar', 31],
                         4 : ['spring', 'apr', 30],
                         5 : ['spring', 'may', 31],
                         6 : ['summer', 'jun', 30],
                         7 : ['summer', 'jul', 31],
                         8 : ['summer', 'aug', 31],
                         9 : ['autumn', 'sep', 30],
                         10: ['autumn', 'oct', 31],
                         11: ['autumn', 'nov', 30]}

    try:
        os.chdir(hysplit_working)

        # Create output directory if it doesn't already exist
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        for year in years:

            isleap = False
            if year % 4 == 0:
                if year % 100 == 0:
                    if year % 400 == 0:
                        isleap = True
                else:
                    isleap = True
            yr = str(year)[-2:]

            for mon in months:

                season = season_month_days[mon][0]
                monname = season_month_days[mon][1]
                mon_len = season_month_days[mon][2]

                if isleap and mon == 2:
                    mon_len = mon_len + 1

                days = range(1, mon_len + 1)

                # Generate a list of the month's meteorology files
                meteofiles = meteofile_lister(meteo_path, meteo_type, mon,
                                              isleap, year)

                # Find total number of meteorology files
                meteofile_count = len(meteofiles)

                for day, hour, alt in itertools.product(days, hours,
                                                        altitudes):

                    # Remove any existing CONTROL or temp files
                    try_to_remove(os.path.join(hysplit_working, filename))
                    try_to_remove(os.path.join(hysplit_working, basename))

                    # Create new control file
                    control = open(os.path.join(hysplit_working, filename),
                                   'w')

                    # Populate trajectory start information
                    ctrltxt = [yr + " {0:02} {1:02} {2:02}\n".format(mon, day,
                                                                     hour),
                               "1\n",
                               "{0!s} {1!s} {2!s}".format(coordinates[0],
                                                          coordinates[1], alt),
                               '.0\n',
                               "{0!s}\n".format(run),
                               "0\n",
                               "10000.0\n",
                               "{0!s}\n".format(meteofile_count)]

                    for meteofile in meteofiles:
                        directory, fname = os.path.split(meteofile)
                        directory = str(directory).replace('\\', '/')
                        ctrltxt.append("{0}/\n".format(directory))
                        ctrltxt.append("{0}\n".format(fname))

                    ctrltxt.append("./\n")
                    ctrltxt.append("{0}\n".format(basename))

                    # Write and flush to disk
                    control.writelines(ctrltxt)
                    control.flush()

                    # Call HYSPLIT to generate trajectory
                    call("C:\\hysplit4\\exec\\hyts_std")

                    # Create descriptive back trajectory filename
                    new_name = (basename + monname + '{:04}'.format(alt) +
                                season + yr +
                                "{0:02}{1:02}{2:02}".format(mon, day, hour))

                    # Generate a forward trajectory, if specified
                    if isbackward and get_forward:
                        forwards_and_backwards(hysplit_working, basename,
                                               filename, output_dir, new_name,
                                               meteofiles)
                    # Check that file of same name isn't in destination already
                    try_to_remove(os.path.join(output_dir, new_name))

                    # Move the trajectory file to output directory
                    os.rename(os.path.join(hysplit_working, basename),
                              os.path.join(output_dir, new_name))

                    if get_clippedtraj:
                        clip_traj(output_dir, new_name)

    finally:
        os.chdir(orig_dir)

    return None


def forwards_and_backwards(hysplit_working, backtraj_fname, control_fname,
                           output_dir, new_name, meteofiles):
    """
    Takes a back trajectory and calculates a forward trajectory from its
        location at the earliest point in time.

    Parameters
    ----------
    hysplit_working : string
        Full or relative path to HYSPLIT working directory
    backtraj_fname : string
        Basename of current back trajectory file
    control_fname : string
        Name of the control file, which should be 'CONTROL'
    output_dir : string
        Full or relative path to back trajectory output directory
    new_name : string
        Back trajectory basename plus month, altitude, season, YYMMDDHH
    meteofiles : list of strings
        List of full paths to meteorology files used to construct
        back trajectory

    Returns
    -------
    Nothing, HYSPLIT is called and results are stored in output_dir/forwardtraj

    """

    # Initialize
    meteofile_count = len(meteofiles)
    last_step_data = []
    forward_fname = new_name + 'FORWARD'

    trajheader = ['Year',
                  'Month',
                  'Date',
                  'Hour (UTC)',
                  'Time step (hr)',
                  'Latitude',
                  'Longitude',
                  'Altitude (magl)']

    # Open the file
    backtraj = open(backtraj_fname, 'r')
    # Get through meteofile info at the head of the file
    while True:
        line = backtraj.readline()
        if 'PRESSURE' in line:
            break

    # Acquire rest of header
    header = line.split()[1:]
    columns = 10 + 2 + len(header)

    # Check if data is multiline
    # Multiline data not supported for some applications (clustering)
    if columns > 20:
        multiline = True
        line1 = backtraj.readline()
        line2 = backtraj.readline()
        linelength = len(line1) + len(line2)
    else:
        multiline = False
        line = backtraj.readline()
        linelength = len(line)

    lastline_start = (linelength + 1) * -1

    # Skip to the last line, read it in
    backtraj.seek(lastline_start, 2)
    last_step = backtraj.readline().split()
    if multiline:
        last_step.extend(backtraj.readline().split())

    # Close file
    backtraj.close()

    # Get forward trajectory start information
    last_step_data.extend(last_step[2:6])
    last_step_data.extend(last_step[8:12])

    year = int(last_step_data[trajheader.index('Year')])
    mon  = int(last_step_data[trajheader.index('Month')])
    day  = int(last_step_data[trajheader.index('Date')])
    hour = int(last_step_data[trajheader.index('Hour (UTC)')])
    lat  = float(last_step_data[trajheader.index('Latitude')])
    lon  = float(last_step_data[trajheader.index('Longitude')])
    alt  = float(last_step_data[trajheader.index('Altitude (magl)')])
    run  = (float(last_step_data[trajheader.index('Time step (hr)')])) * -1

    # Sometimes startheight is greater than modeltop at 10000 m
    if alt >= 10000:
        alt = 9999

    output_fdir = os.path.join(output_dir, 'forwardtraj')

    # Make forward trajectory repository if it doesn't exist
    if not os.path.isdir(output_fdir):
        os.mkdir(output_fdir)

    # Remove (if present) any existing CONTROL or temp files
    try_to_remove(os.path.join(hysplit_working, control_fname))
    try_to_remove(os.path.join(hysplit_working, forward_fname))
    try_to_remove(os.path.join(output_fdir, forward_fname))

    # Make a new control file and a new temp file
    control = open(os.path.join(hysplit_working, control_fname), 'w')

    # Populate what we want in the file
    ctrltxt = ["{0:02} {1:02} {2:02} {3:02}\n".format(year, mon, day, hour),
               "1\n",
               "{0!s} {1!s} {2!s}\n".format(lat, lon, alt),
               "{0!s}\n".format(run),
               "0\n",
               "10000.0\n",
               "{0!s}\n".format(meteofile_count)]

    for meteorology_file in meteofiles:
        directory, fname = os.path.split(meteorology_file)
        directory = str(directory).replace('\\', '/')
        ctrltxt.append("{0}/\n".format(directory))
        ctrltxt.append("{0}\n".format(fname))

    ctrltxt.append("./\n")
    ctrltxt.append("{0}\n".format(forward_fname))

    # Write and flush to disk
    control.writelines(ctrltxt)
    control.flush()

    # Call HYSPLIT to generate the next trajectory
    call("C:\\hysplit4\\exec\\hyts_std")

    # Move the trajectory file to the desired output directory
    os.rename(os.path.join(hysplit_working, forward_fname),
              os.path.join(output_fdir, forward_fname))

    return None


def clip_traj(output_dir, new_name):
    """
    Creates a new back trajectory file with all meteorological data
        clipped off, perfect for clustering multiline files.

    New file lives in subdirectory in output_dir

    Parameters
    ----------
    output_dir : string
        Full or relative path to back trajectory output directory
    new_name : string
        Back trajectory basename plus month, altitude, season, YYMMDDHH

    Returns
    -------
    None, creates a file

    """

    # Open back trajectory file
    originaltraj = open(os.path.join(output_dir, new_name), 'r')

    # Initialize
    part1 = []
    part2 = []
    clipped_fname = new_name + 'CLIPPED'

    # Get file header information
    while True:
        line = originaltraj.readline()

        if 'PRESSURE' in line:
            break

        part1.append(line)

    # See if data is multiline or not
    datacolumns = line.split()[1:]
    num_datacolumns = 12 + len(datacolumns)

    multiline = False
    if num_datacolumns > 20:
        multiline = True

    # Put in data header line
    line = '     1 PRESSURE\n'
    part1.append(line)

    while True:
        line = originaltraj.readline()

        if line == '':
            part2.append(line)
            break

        if multiline:
            line = line + originaltraj.readline()

        # Put clipped data lines into list
        part2.append(line[:92] + '\n')

    originaltraj.close()

    output_cdir = os.path.join(output_dir, 'clippedtraj')

    if not os.path.isdir(output_cdir):
        os.mkdir(output_cdir)

    # Remove file if it exists
    try_to_remove(os.path.join(output_cdir, clipped_fname))

    clippedtraj = open(os.path.join(output_cdir, clipped_fname), 'w')

    clippedtraj.writelines(part1)
    clippedtraj.writelines(part2)

    clippedtraj.flush()
    clippedtraj.close()

    return None


def try_to_remove(string):
    # Will error if no such file exists. We just want to make sure there is
    # no file there, so try first - if it fails, continue as normal.
    try:
        os.remove(string)
    except OSError:
        pass


def meteofile_lister(meteo_path, meteo_type, mon, isleap, year):
    """
    Takes a month and information about leap year status to create a list of
        meteorology files to be output to CONTROL.

    Parameters
    ----------
    meteo_path : string
        first part of meteorology file path name
        ex. 'C:/hysplit4/working/gdas1.'
    meteo_type : string
        ['gdas1'|'era-interim'].  Type of ARL-formatted data files provided.
    mon : int
        int representing the month
    isleap : Boolean
        [True|False].  True if year of simulation is a leap year
    year : int
        int representing the year
    backward : Boolean
        [True|False].  True if performing backwards trajectory calculation

    Returns
    -------
    meteofiles : list of strings
        List of strings representing complete file paths
    """

    fpath = meteo_path + '/'

    if meteo_type == 'gdas1':

        if year < 2005:
            raise ValueError('Year out of range of GDAS1 dataset!')

        f = 'gdas1.'

        # Initialize lists of parts of filename strings
        months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun',
                  'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
        weeks = ['.w1', '.w2', '.w3', '.w4', '.w5']

        # Initialize empty lists
        leadweek = []
        midweeks = []
        endweeks = []
        meteofiles = []

        # Assemble filenames for first and second week of next month
        # If it is December, weeks are from January of next year
        if (mon == 12):
            for week in weeks[:2]:
                fname = f + months[0] + str(year + 1)[-2:] + week
                leadweek.append(fpath + fname)

        # Otherwise, from the next month
        else:
            for week in weeks[:2]:
                fname = f + months[mon] + str(year)[-2:] + week
                leadweek.append(fpath + fname)

        # Assemble filenames for weeks of this month
        # If it is February and not a leap year, there are only four week files
        if (mon == 2) & (isleap is False):
            for week in weeks[0:4]:
                fname = f + months[mon - 1] + str(year)[-2:] + week
                midweeks.append(fpath + fname)

        # Otherwise, there are five week files
        else:
            for week in weeks:
                fname = f + months[mon - 1] + str(year)[-2:] + week
                midweeks.append(fpath + fname)

        # Assemble filename(s) for last week(s) of previous month
        # If it is March and it is not a leap year, grab only the 4th week file
        if (mon == 3) & (isleap is False):
            fname = f + months[mon - 2] + str(year)[-2:] + weeks[3]
            endweeks.append(fpath + fname)

        # If it is January, need the 4-5th week files from previous year's Dec
        elif mon == 1:
            for week in weeks[3:]:
                fname = f + months[11] + str(year - 1)[-2:] + week
                endweeks.append(fpath + fname)

        # Otherwise, get the 4-5th week files from previous month
        else:
            for week in weeks[3:]:
                fname = f + months[mon - 2] + str(year)[-2:] + week
                endweeks.append(fpath + fname)

        meteofiles.extend(leadweek)
        meteofiles.extend(midweeks)
        meteofiles.extend(endweeks)

        for meteofile in meteofiles:
            if not os.path.exists(meteofile):
                raise OSError('Meteorology file does not exist!')

    elif meteo_type is 'era-interim':

        # Check that meteorology exists for time range
        if year < 1979:
            raise ValueError('Year out of range of ERA-interim dataset!')

        f = 'ERA'
        meteofiles = []

        # Get files from the actual month
        p1 = fpath + f + str(year) + '_' + '{:02}'.format(mon) + '_1'
        p2 = fpath + f + str(year) + '_' + '{:02}'.format(mon) + '_2'
        p3 = fpath + f + str(year) + '_' + '{:02}'.format(mon) + '_3'

        # Get files from months before and after
        # Special cases January and December
        if mon == 1:
            before = fpath + f + str(year - 1) + '_12_3'
        else:
            before = (fpath + f + str(year) + '_' +
                      '{:02}'.format(mon - 1) + '_3')

        if mon == 12:
            after = fpath + f + str(year + 1) + '_01_1'
        else:
            after = (fpath + f + str(year) + '_' +
                     '{:02}'.format(mon + 1) + '_1')

        # Put into list
        meteofiles.append(before)
        meteofiles.append(p1)
        meteofiles.append(p2)
        meteofiles.append(p3)
        meteofiles.append(after)

        # Check that meteorology exists
        for meteofile in meteofiles:
            if not os.path.exists(meteofile):
                raise OSError('Meteorology file does not exist!')

    return meteofiles
