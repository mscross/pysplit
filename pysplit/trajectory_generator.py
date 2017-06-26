from __future__ import division
import os
from subprocess import call
import itertools
import fnmatch
from calendar import monthrange


def generate_bulktraj(basename, hysplit_working, output_dir, meteo_dir, years,
                      months, hours, altitudes, coordinates, run,
                      monthslice=slice(0, 32, 1), meteo_bookends=([4, 5], [1]),
                      get_reverse=False, get_clipped=False,
                      hysplit="C:\\hysplit4\\exec\\hyts_std"):
    """
    Generate sequence of trajectories within given time frame(s).

    Run bulk sequence of HYSPLIT simulations over a given time and at different
    altitudes (likely in meters above ground level).  Uses either weekly or
    semi-monthly data with the filename format of *mon*YY*#.
    Results are written to ``output_dir``.

    This does not set along-trajectory meteorological output- edit SETUP.CFG
    in the HYSPLIT working directory or in the HYSPLIT4 GUI to reflect
    desired output variables.  It is also recommended to change TRATIO
    in the SETUP file to 0.25 to limit integration error.

    Parameters
    ----------
    basename : string
        Base for all files output in this run
    hysplit_working : string
        Full or relative path to the HYSPLIT working directory.
    output_dir : string
        Full or relative path to the desired output directory.
    meteo_dir : string
        Full or relative path to the location of the meteorology files.
    years : list of ints
        The year(s) to run simulations
    months : list of ints
        The month(s) to run simulations
    hours : list of ints
        Parcel launching times in UTC.
    altitudes : list of ints
        The altitudes (usually meters above ground level) from which
        parcels will be launched.  Must be less than model top (10000 m)
    coordinates : tuple of floats
        The parcel (latitude, longitude) launch location in decimal degrees.
    run : int
        Length in hours of simulation.  To calculate back trajectories,
        ``run`` must be negative.
    monthslice : slice object
        Default slice(0, 32, 1).  Slice to apply to range of days in month.
        Use to target particular day or range of days, every x number of days,
        etc.  NOTE: slice is 0 indexed, days start with 1.
    meteo_bookends : tuple of lists of ints
        Default ([4, 5], [1]).  To calculate a month of trajectories, files
        from the previous and month must be included.  The default is optimized
        for weekly meteorology and indicates that weeks 4 and 5 from the
        previous month and the first week of the next month must be included
        to run the entire current month of trajectories.  The user is
        responsible for making sure the correct bookends for their trajectory
        length and meteorology file periods are provided.
    get_reverse : Boolean
        Default ``False``.  If ``True``, then from the last point of each
        trajectory a new parcel will be launched in the opposite direction.
        These reverse trajectories are stored in a subfolder in ``output_dir``
    get_clipped : Boolean
        Default ``False``.   If ``True``, takes a trajectory file and
        outputs a version of the file containing only path information.
        Provided to support clustering of trajectories with multiline data,
        which was produced in HSYPLI versions prior to January 2017 (854)
        when more than 7 along-trajectory output variables were selected.
    hysplit : string
        Default "C:\\hysplit4\\exec\\hyts_std".  The location of the "hyts_std"
        executable that generates trajectories.  This is the default location
        for a typical PC installation of HYSPLIT

    """
    controlfname = 'CONTROL'

    # Get directory information, make directories if necessary
    cwd = os.getcwd()

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    output_rdir = os.path.join(output_dir, 'reversetraj')
    output_cdir = os.path.join(output_dir, 'clippedtraj')
    meteo_dir = meteo_dir.replace('\\', '/')

    if get_reverse and not os.path.isdir(output_rdir):
        os.mkdir(os.path.join(output_rdir))

    if get_clipped and not os.path.isdir(output_cdir):
        os.mkdir(os.path.join(output_cdir))

    # Initialize dictionary of months, seasons
    n_hemisphere = True
    if coordinates[0] < 0:
        n_hemisphere = False

    mon_dict = _mondict(n_hem=n_hemisphere)

    try:
        os.chdir(hysplit_working)

        # Iterate over years and months
        for y, m in itertools.product(years, months):
            season = mon_dict[m][0]
            m_str = mon_dict[m][1]
            m_len = monthrange(y, m)[1]

            days = range(1, m_len + 1)[monthslice]

            # Assemble list of meteorology files
            meteofiles = _meteofinder(meteo_dir, meteo_bookends, m, y,
                                      mon_dict)

            yr = _year2string(y)

            # Iterate over days, hours, altitudes
            for d, h, a in itertools.product(days, hours, altitudes):

                # Add timing and altitude to basename to create unique name
                trajname = (basename + m_str + '{:04}'.format(a) + season +
                            yr + "{0:02}{1:02}{2:02}".format(m, d, h))

                final_trajpath = os.path.join(output_dir, trajname)

                # Remove any existing CONTROL or temp files
                _try_to_remove(controlfname)
                _try_to_remove(trajname)
                _try_to_remove(final_trajpath)

                # Populate CONTROL file with trajectory initialization data
                _populate_control(coordinates, yr, m, d, h, a, meteo_dir,
                                  meteofiles, run, controlfname, trajname)

                # Call executable to calculate trajectory
                call(hysplit)

                # Generate reverse and/or clipped trajectories, if indicated
                if get_reverse:
                    _reversetraj_whilegen(trajname, run, hysplit, output_rdir,
                                          meteo_dir, meteofiles, controlfname)

                if get_clipped:
                    _cliptraj(output_cdir, trajname)

                # Move the trajectory file to output directory
                os.rename(trajname, final_trajpath)

    # Revert current working directory
    finally:
        os.chdir(cwd)

def generate_traj(basename, hysplit_working, output_dir, meteo_dir, years,
                          months, hours, altitudes, coordinates, run,
                          monthslice=slice(0, 32, 1), meteo_bookends=([4, 5], [1]),
                          get_reverse=False, get_clipped=False,
                          hysplit="C:\\hysplit4\\exec\\hyts_std"):

    """
    Generate trajectory with given time frame.
    
    Run a HYSPLIT simulation over a given time and 
    altitude (likely in meters above ground level).  Uses either weekly or
    semi-monthly data with the filename format of *mon*YY*#.
    Results are written to ``output_dir``.

    This does not set along-trajectory meteorological output- edit SETUP.CFG
    in the HYSPLIT working directory or in the HYSPLIT4 GUI to reflect
    desired output variables.  It is also recommended to change TRATIO
    in the SETUP file to 0.25 to limit integration error.

    Parameters
    ----------
    basename : string
        Base for all files output in this run
    hysplit_working : string
        Full or relative path to the HYSPLIT working directory.
    output_dir : string
        Full or relative path to the desired output directory.
    meteo_dir : string
        Full or relative path to the location of the meteorology files.
    year : 
        The year to run simulation
    month : 
        The month to run simulation
    hour : 
        Parcel launching time in UTC.
    altitude :
        The altitude (usually meters above ground level) from which
        parcel will be launched.  Must be less than model top (10000 m)
    coordinates : tuple of floats
        The parcel (latitude, longitude) launch location in decimal degrees.
    run : int
        Length in hours of simulation.  To calculate back trajectories,
        ``run`` must be negative.
    monthslice : slice object
        Default slice(0, 32, 1).  Slice to apply to day in month.
        Use to target particular day,
        etc.  NOTE: slice is 0 indexed, days start with 1.
    meteo_bookends : tuple of lists of ints
        Default ([4, 5], [1]).  To calculate a trajectory, files
        from the previous and month must be included.  The default is optimized
        for weekly meteorology and indicates that weeks 4 and 5 from the
        previous month and the first week of the next month must be included
        to run the entire current month of trajectories.  The user is
        responsible for making sure the correct bookends for their trajectory
        length and meteorology file periods are provided.
    get_reverse : Boolean
        Default ``False``.  If ``True``, then from the last point of
        trajectory a new parcel will be launched in the opposite direction.
        These reverse trajectories are stored in a subfolder in ``output_dir``
    get_clipped : Boolean
        Default ``False``.   If ``True``, takes a trajectory file and
        outputs a version of the file containing only path information.
        Provided to support clustering of trajectories with multiline data,
        which was produced in HSYPLI versions prior to January 2017 (854)
        when more than 7 along-trajectory output variables were selected.
    hysplit : string
        Default "C:\\hysplit4\\exec\\hyts_std".  The location of the "hyts_std"
        executable that generates trajectories.  This is the default location
        for a typical PC installation of HYSPLIT
        
    """

    controlfname = 'CONTROL'
    # Get directory information, make directories if necessary
    cwd = os.getcwd()

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    output_rdir = os.path.join(output_dir, 'reversetraj')
    output_cdir = os.path.join(output_dir, 'clippedtraj')
    meteo_dir = meteo_dir.replace('\\', '/')

    if get_reverse and not os.path.isdir(output_rdir):
        os.mkdir(os.path.join(output_rdir))

    if get_clipped and not os.path.isdir(output_cdir):
        os.mkdir(os.path.join(output_cdir))

    # Initialize dictionary of months, seasons

    n_hemisphere = True
    if coordinates[0] < 0:
        n_hemisphere = False

    mon_dict = _mondict(n_hem=n_hemisphere)

    try:
        os.chdir(hysplit_working)

        # Iterate over years and months

        for y, m in zip(years, months):
            season = mon_dict[m][0]
            m_str = mon_dict[m][1]
            m_len = monthrange(y, m)[1]
            
            days = range(1, m_len + 1)[monthslice]

            # Assemble list of meteorology files
            meteofiles = _meteofinder(meteo_dir, meteo_bookends, m, y,
                                      mon_dict)

            yr = _year2string(y)

            # Iterate over days, hours, altitudes
            for d, h, a in zip(days, hours, altitudes):
              
                # Add timing and altitude to basename to create unique name
                trajname = (basename + m_str + '{:04}'.format(a) + season +
                            yr + "{0:02}{1:02}{2:02}".format(m, d, h))

                final_trajpath = os.path.join(output_dir, trajname)
                
                # Remove any existing CONTROL or temp files

                _try_to_remove(controlfname)
                _try_to_remove(trajname)
                _try_to_remove(final_trajpath)

                # Populate CONTROL file with trajectory initialization data
                _populate_control(coordinates, yr, m, d, h, a, meteo_dir,
                                  meteofiles, run, controlfname, trajname)

                # Call executable to calculate trajectory
                call(hysplit)

                # Generate reverse and/or clipped trajectories, if indicated
                if get_reverse:
                    _reversetraj_whilegen(trajname, run, hysplit, output_rdir,
                                          meteo_dir, meteofiles, controlfname)

                if get_clipped:
                    _cliptraj(output_cdir, trajname)

                # Move the trajectory file to output directory
                os.rename(trajname, final_trajpath)

    # Revert current working directory
    finally:
        os.chdir(cwd)
        
        
def _reversetraj_whilegen(trajname, run, hysplit, output_rdir, meteo_dir,
                          meteofiles, controlfname):
    """
    Calculate reverse trajectory during main trajectory generation.

    Calculates a new trajectory ('reverse trajectory') from the endpoint of the
    trajectory just calculated in the main generation sequence and running
    in the opposite direction.

    Parameters
    ----------
    trajname : string
        The file name of the just-calculated trajectory.  New backwards
        trajectory will be named ``trajname`` + 'REVERSE'
    run : int
        The length in hours of the trajectory simulation
    hysplit : string
        The location of the executable that calculates trajectories
    output_rdir : string
        The subdirectory for reverse trajectories.
    meteo_dir : string
        The location of the meteorology files
    meteofiles : string
        The list of meteorology files required to calculate the
        reverse trajectory.
    controlfname : string
        The name of the control file, which should be 'CONTROL'

    """
    # Initialize name, path
    reversetrajname = trajname + 'REVERSE'
    final_rtrajpath = os.path.join(output_rdir, reversetrajname)

    with open(trajname) as traj:
        contents = traj.readlines()

        last_timepoint = -1

        # for multiline files, critical information is in contents[-2]
        if len(contents[-1]) < len(contents[-2]):
            last_timepoint = -2

        data = contents[last_timepoint].split()

    # Get reverse trajectory start information
    year = int(data[2])
    mon  = int(data[3])
    day  = int(data[4])
    hour = int(data[5])
    lat  = float(data[9])
    lon  = float(data[10])
    alt  = float(data[11])
    run  = run * -1

    # Sometimes start height is greater than 10000 m model top
    if alt >= 10000:
        alt = 9999

    yr = '{:02}'.format(year)

    # Remove (if present) any existing CONTROL or temp files
    _try_to_remove(controlfname)
    _try_to_remove(reversetrajname)
    _try_to_remove(final_rtrajpath)

    # Populate control text
    _populate_control((lat, lon), yr, mon, day, hour, alt, meteo_dir,
                      meteofiles, run, controlfname, reversetrajname)

    # Call executable
    call(hysplit)

    # Move the trajectory file to the desired output directory
    os.rename(reversetrajname, final_rtrajpath)


def _cliptraj(output_cdir, trajname):
    """
    Create clipped trajectory file from original file.

    Creates a new trajectory file containing only header and path information
    from a newly generated trajectory.  Only necessary if files are multiline.

    Parameters
    ----------
    output_cdir : string
        Full or relative path to clipped trajectory output directory
    trajname : string
        Name of trajectory file to clip.  New file will be named ``trajname`` +
        'CLIPPED'

    """
    # Initialize name, path, data list
    clippedtrajname = trajname + 'CLIPPED'
    final_ctrajpath = os.path.join(output_cdir, clippedtrajname)

    clipdata = []

    with open(trajname) as original:

        # Read in all lines of file
        contents = original.readlines()

        # Initialize markers
        skip = False
        atdata = False
        multiline = False

        # Iterate through lines
        for ind, line in enumerate(contents):
            # Skip line only triggered if multiline, after at data
            if skip:
                skip = False
                continue

            # Once at data, only need first 92 char of line(s), append data
            if atdata:
                clipdata.append(line[:92] + '\n')
                if multiline:
                    skip = True
                continue

            # PRESSURE marker tripped first
            if 'PRESSURE' in line:
                if len(contents[ind + 1]) > len(contents[ind + 2]):
                    multiline = True
                # Append last line of header, now at data
                clipdata.append(line[:15] + '\n')
                atdata = True
                continue

            # Append header data as is
            clipdata.append(line)

    # Get rid of temporary files and files with the same path
    _try_to_remove(clippedtrajname)
    _try_to_remove(final_ctrajpath)

    # Write data to file and move
    with open(clippedtrajname, 'w') as ctraj:
        ctraj.writelines(clipdata)

    os.rename(clippedtrajname, final_ctrajpath)


def _meteofinder(meteo_dir, meteo_bookends, mon, year, mon_dict):
    """
    Get list of meteorology files.

    Creates list of files in storage location ``meteo_dir`` that belong
    to the given month and year, plus the necessary files from previous
    and the next months (``meteo_bookends``).

    For successful meteofinding, separate different meteorology types into
    different folders and name weekly or semi-monthly files according to the
    following convention:
        *mon*YY*#
    where the * represents a Bash wildcard.

    Parameters
    ----------
    meteo_dir : string
        Full or relative path to the location of the meteorology files.
    meteo_bookends : tuple of lists of ints
        To calculate a month of trajectories, files from the previous and next
        month must be included.  This indicates which file numbers from the
        previous month and which from the next month are necessary.
        The user is responsible for making sure the correct bookends for their
        trajectory length and meteorology file periods are provided.
    mon : int
        The integer representation of the current month.  Converted to a
        3-letter string to find meteorology files.
    year : int
        The integer representation of the current year.  Converted to a length
        2 string to find meteorology files.
    mon_dict : dictionary
        Dictionary keyed by month integer, with lists of [season, mon]

    Returns
    -------
    meteofiles : list of strings
        List of strings representing the names of the required
        meteorology files

    """
    # Current working directory set in generate_bulktraj() environment
    orig_dir = os.getcwd()

    # Initialize lists, count
    meteofiles = []
    file_number = -1

    # Get the strings that will match files for the previous, next,
    # and current months
    prv, nxt, now = _monyearstrings(mon, year, mon_dict)

    # Change directory and walk through files
    try:
        os.chdir(meteo_dir)

        _, _, files = next(os.walk('.'))

        # Order of files to CONTROL doesn't matter
        for each_file in files:
            if fnmatch.fnmatch(each_file, now):
                meteofiles.append(each_file)
            elif fnmatch.fnmatch(each_file, prv):
                if int(each_file[file_number]) in meteo_bookends[0]:
                    meteofiles.append(each_file)
            elif fnmatch.fnmatch(each_file, nxt):
                if int(each_file[file_number]) in meteo_bookends[1]:
                    meteofiles.append(each_file)

    finally:
        os.chdir(orig_dir)

    if len(meteofiles) == 0:
        raise OSError('0 files found for month/year %(mon)d / %(year)d'
                      %{'mon': mon, 'year': year})

    return meteofiles


def _populate_control(coords, year, month, day, hour, alt,
                      meteo_dir, meteofiles, run, controlfname, trajname):
    """
    Initialize and write CONTROL text to file (called CONTROL).

    Parameters
    ----------
    coordinates : tuple of floats
        The parcel (latitude, longitude) launch location in decimal degrees.
    years : list of ints
        The year of the simulation
    months : list of ints
        The month of the simulation
    hours : list of ints
        Parcel launching times in UTC.
    alt : int
        The altitude (usually meters above ground level) from which
        parcel will be launched.  Must be less than model top (10000 m)
    meteo_dir : string
        Full or relative path to the location of the meteorology files.
    meteofiles : list of strings
        List of strings representing the names of the required
        meteorology files
    run : int
        Length in hours of simulation.
    controlfname : string
        The name of the control file, which should be 'CONTROL'
    trajname : string
        The intended name of the trajectory file

    """
    controltext = [year + " {0:02} {1:02} {2:02}\n".format(month, day, hour),
                   "1\n",
                   "{0!s} {1!s} {2!s}\n".format(coords[0], coords[1], alt),
                   "{0!s}\n".format(run),
                   "0\n",
                   "10000.0\n",
                   "{0!s}\n".format(len(meteofiles))]

    for fname in meteofiles:
        controltext.append("{0}/\n".format(meteo_dir))
        controltext.append("{0}\n".format(fname))

    controltext.append("./\n")

    controltext.append("{0}\n".format(trajname))

    with open(controlfname, 'w') as control:
        control.writelines(controltext)


def _year2string(year):
    """
    Helper function, takes a four digit integer year, makes a length-2 string.

    Parameters
    ----------
    year : int
        The year.

    Returns
    -------
    Length-2 string representation of ``year``

    """
    return '{0:02}'.format(year % 100)


def _monyearstrings(mon, year, mon_dict):
    """
    Increment the months and potentially the years.

    Assemble the strings that will allow ``_meteofinder`` to get correct files.

    Parameters
    ----------
    mon : int
        Integer representation of the month
    year : int
        Integer representation of the year
    mon_dict : dictionary
        Dictionary keyed by month integer, with lists of [season, mon]

    Returns
    -------
    prv : string
        Signature for gathering the meteorology files for the previous month
    nxt : string
        Signature for gathering the meteorology files for the next month
    now : string
        Signature for gathering the meteorology files for the current month

    """
    next_year = year
    prev_year = year

    next_mon = mon + 1
    prev_mon = mon - 1

    if prev_mon == 0:
        prev_mon = 12
        prev_year = year - 1
    if next_mon == 13:
        next_mon = 1
        next_year = year + 1

    w = '*'

    prv = w + mon_dict[prev_mon][1] + w + str(prev_year) + w
    nxt = w + mon_dict[next_mon][1] + w + str(next_year) + w
    now = w + mon_dict[mon][1] + w + str(year) + w

    return prv, nxt, now


def _mondict(n_hem=True):
    """
    Get a dictionary of season and month string.

    Parameters
    ----------
    n_hem : Boolean
        Default True.  Indicates hemisphere of parcel launch and thus
        actual season.

    Returns
    -------
    season_month_dict : dictionary
        Dictionary keyed by month integer, with lists of [season, mon]

    """
    if n_hem:
        season_month_dict = {12: ['winter', 'dec'],
                             1 : ['winter', 'jan'],
                             2 : ['winter', 'feb'],
                             3 : ['spring', 'mar'],
                             4 : ['spring', 'apr'],
                             5 : ['spring', 'may'],
                             6 : ['summer', 'jun'],
                             7 : ['summer', 'jul'],
                             8 : ['summer', 'aug'],
                             9 : ['autumn', 'sep'],
                             10: ['autumn', 'oct'],
                             11: ['autumn', 'nov']}
    else:
        season_month_dict = {12: ['summer', 'dec'],
                             1 : ['summer', 'jan'],
                             2 : ['summer', 'feb'],
                             3 : ['autumn', 'mar'],
                             4 : ['autumn', 'apr'],
                             5 : ['autumn', 'may'],
                             6 : ['winter', 'jun'],
                             7 : ['winter', 'jul'],
                             8 : ['winter', 'aug'],
                             9 : ['spring', 'sep'],
                             10: ['spring', 'oct'],
                             11: ['spring', 'nov']}

    return season_month_dict


def _try_to_remove(string):
    """
    Check if file exists, and either remove it or pass.

    Parameters
    ----------
    string : string
        Name of file to attempt to remove

    """
    try:
        os.remove(string)
    except OSError:
        pass


def _day2filenum(interval, day):
    """
    Convert a date to corresponding file number.

    Results depend on file interval- weekly, daily, semi-monthly.

    Parameters
    ----------
    interval : string
        The file interval.  Daily, weekly, or semi-monthly accepted,
        represented by lower case first letter.
    day : string
        A number indicating the date.

    Returns
    -------
    filenum : string
        The number of the file within the month of meteorology.
    """
    if interval == 'w':
        filenum = str(((int(day) - 1) // 7) + 1)
    elif interval == 's':
        filenum = str(((int(day) - 1) // 15) + 1)
    elif interval == 'd' or interval == 'm':
        filenum = day
    else:
        raise ValueError('Meteorology interval not recognized')

    return filenum
