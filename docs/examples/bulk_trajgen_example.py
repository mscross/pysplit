"""
=====================
Trajectory Generation
=====================

PySPLIT includes a method that calls HYSPLIT to generate trajectories
launched from a single location over as many months and years as desired in
a single comprehensive batch.


Requirements
------------

HYSPLIT installation
~~~~~~~~~~~~~~~~~~~~

Unlike the majority of PySPLIT, this feature requires a local installation
of HYSPLIT, which can be acquired from the NOAA ARL
here<http://ready.arl.noaa.gov/HYSPLIT.php>_.

The location of the `hyts_std` executable must be known.

Meteorology
~~~~~~~~~~~

Meteorology data sets corresponding to the desired time and geospatial
extent must be available to successfully generate trajectories.
These data sets must be in the ARL-packed format.  All meteorology downloaded
through the HYSPLIT ARL data FTP is already in this format, and HYSPLIT
contains some utilities to convert other file types to this format. 
Acquistion of appropriate files and the nature of
the ARL-packed format and conversion is beyond the scope of this example,
and for more information consult the official `User's
Guide<http://www.arl.noaa.gov/documents/reports/hysplit_user_guide.pdf`_
if necessary before proceeding.

Prior to attempting trajectory generation with PySPLIT, the user should ensure
that meteorology file names follow a format of '*mon*YY*#*', where:

* '*' is a Bash-style wildcard
* 'mon' is a three-letter lower-case abbreviation of the month
* 'YY' is a two digit integer representation of the year
* '#' is the number of the file within the month

For example, '*jan*07*2' will match files named 'gdas1.jan2007.w2' and
'edas.jan2007_2'.

This naming convention is required because PySPLIT examines file names to
determine which meteorology files correspond to the desired date of
trajectory launch.  It is strongly recommended that users keep meteorology
files of different types/origins, like the two files in the example above,
in separate directories as PySPLIT does not differentiate between the two.

Output along-trajectory meteorology variables may be selected by interacting
with HYSPLIT.


Generating Trajectories
-----------------------

As always, we begin by importing the package.
"""

import pysplit

"""
For clarity, we will define all required arguments outside of the method.

The first three arguments indicate locations of the HYSPLIT working directory,
the desired and existing trajectory storage directory, and the directory
that contains the meteorology files.  This example uses 1-degree weekly
GDAS archive meteorology downloaded via FTP.
"""

working_dir = r'C:/hysplit4/working'
storage_dir = r'C:/trajectories/colgate'
meteo_dir = r'E:/gdas'

"""
The next argument is the basename of each trajectory.  Each trajectory
filename will consist of this basename followed by the altitude and season,
and then the year, month, day, and hour in the format YYMMDDHH.
"""

basename = 'colgate'

"""
The following arguments are the initialization information for the 
trajectory generation run:  lists or ranges of years, months, hours (in UTC,
NOT local time), and altitudes in meters above ground level.
We also supply the tuple of initialization coordinates in decimal degrees,
which here indicate Colgate University in Hamilton, NY, and the integer run
time of the trajectories.  The trajectory direction in time is specified by
the sign of the run time (negative is backwards).
"""

years = [2007, 2011]
months = [1, 8]
hours = [11, 17, 23]
altitudes = [500, 1000, 1500]
location = (42.82, -75.54)
runtime = -120

"""
There are five keyword arguments in this method.  The first is a slice object
that is applied to the range of days in each month, allowing us to generate
trajectories for all (default ``slice(0, 32, 1)) or a subset of days in each
month.  Here we choose to launch trajectories every other day. 

We have so far set up our run such that trajectories are launched every other
day in January and August, and their paths are calculated thence backwards
five days.  For trajectories launched at the beginning of a month, therefore,
we need meteorology data from the previous month.  The ``meteo_bookends``
indicates which meteorology files from the previous and next months are
needed.  As we are using weekly datasets, we need
weeks 4 (days 22-28) and weeks 5.  The default value (``[[4, 5], [1]]``)
also includes the first week of next month, but that's ok to include.

The keywords ``get_reverse`` (default False) and ``get_clipped`` (default
True) control whether additional trajectory files are created.  A 'reverse'
trajectory is one launched from the endpoint of an original trajectory, and
run in the opposite direction.  A 'clipped' trajectory does not require
additional calculation, but is a copy of the original trajectory file
containing only path information.  We will acquire both reverse and clipped
trajectories.

The final keyword argument is the location of the HYSPLIT trajectory
executable.  On Windows systems, this will usually be
'C:\\hysplit4\\exec\\hyts_std' (default value).  PySPLIT calls this
executable to calculate trajectories.

Let's call our method.  Only keyword arguments with values we have decided 
to change from default are included below.  This call may take several
minutes to complete, depending on your system parameters.
"""

pysplit.generate_bulktraj(basename, working_dir, storage_dir, meteo_dir,
                          years, months, hours, altitudes, location, runtime,
                          monthslice=slice(0, 32, 2), get_reverse=True)

"""
When complete, ``storage_dir``, will contain 578 trajectory files, as well as two
folders, 'reversetraj' and 'clippedtraj', each containing 578 files.
"""
