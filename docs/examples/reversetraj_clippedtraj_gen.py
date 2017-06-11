"""
=====================================================================
Generating Reverse and Clipped Trajectories After Original Experiment
=====================================================================

In the example ``bulk_trajgen_example.py``, the basic PySPLIT procedure
for generating HYSPLIT trajectories was presented.  While generating
bulk trajectories (``pysplit.generate_bulktraj()``), the user may choose
to also output "clipped" trajectory files, which are copies of the original
trajectory files that contain only path information, and/or "reverse"
trajectory files, which are newly calculated trajectories launched from
the endpoint of the original trajectories in the opposite direction in time.

For large experiments, creating these extra files, particularly the
"reverse" trajectories, can significantly inflate the run time.  So,
instead of generating the complementary "reverse" and "clipped"
trajectories during the main trajectory generation, the user may choose
to generate these extras later and with a subset of trajectories.


Setup
-----

Below, we'll generate back trajectories from the University of Minnesota,
forgoing the "reverse" and "clipped" trajectories (``get_reverse=False``,
``get_clipped=False``).  For a detailed explanation of this call, see
``bulk_trajgen_example.py``.  Then, initialize a TrajectoryGroup.

"""
from __future__ import print_function
import pysplit

pysplit.generate_bulktraj('umn', r'C:/hysplit4/working',
                          r'C:/trajectories/umn_example', r'E:/gdas', [2013],
                          range(1,4), [6, 15, 18, 21], [1000], (44.97, -93.23),
                          -96, get_clipped=False)

trajgroup = pysplit.make_trajectorygroup(r'C:/trajectories/umn_example/*')

"""
In the course of our analysis, we might decide that a subset of trajectories
is more appropriate for our purpose.  Let's create a new ``TrajectoryGroup``
with all the trajectories with a temperature at t=0 of 0 degrees C and greater.

"""
warm_trajlist = [traj for traj in trajgroup if traj.data.Temperature_C[0] > 0.0]

warm_trajgroup = pysplit.TrajectoryGroup(warm_trajlist)

"""
The new ``TrajectoryGroup`` is much smaller:

"""
print(trajgroup.trajcount)
print(warm_trajgroup.trajcount)

"""
Generation
----------

Now we generate the reverse trajectories and create the clipped trajectory 
files.

In ``Trajectory.generate_reversetraj()``, the arguments are
the HYSPLIT working directory, and the meteorology file location.  The kwargs
indicate the interval between meteorology files ('monthly', 'semimonthly',
'weekly', 'daily'), the directory to store the reverse trajectories in, if 
some location other than the default (in a folder in the trajectory directory)
is desired; and the location of the hysplit executable.

In ``Trajectory.generate_clippedtraj()``, the kwarg is the storage
location for the clipped trajectory files, if some location other than the default
is desired.  We use the default locations for clipped and reverse trajectories
in this example.

If any clipped or reverse trajectories corresponding to the trajectories
in ``warm_trajgroup`` exist in the storage directories, they will be overwritten.

If the reverse and clipped directories don't exist, they will be
created. 

"""
for traj in warm_trajgroup:
    traj.generate_reversetraj(r'C:/hysplit4/working', r'E:/gdas',
                              meteo_interval='weekly',
                              hysplit="C:\\hysplit4\\exec\\hyts_std")

    traj.generate_clippedtraj()

"""
Loading and Mapping Reverse Trajectories
----------------------------------------
Load the reverse trajectory data via ``Trajectory.load_reversetraj()``.

For details on the sample basic plotting procedure, see
``basic_plotting_example.py``.  The original trajectories are plotted in
gold and the reverse trajectories in maroon.  Ideally, the reverse
trajectories will overlay the original trajectories.

"""
for traj in warm_trajgroup:
    traj.load_reversetraj()
    
mapcorners =  [-150, 15, -50, 65]
standard_pm = None

testmap = pysplit.MapDesign(mapcorners, standard_pm)

bmap = testmap.make_basemap()

for traj in warm_trajgroup[::5]:
    bmap.plot(*traj.path.xy, c='#FFB71E', latlon=True, zorder=20)
    bmap.plot(*traj.path_r.xy, c='#5B0013', latlon=True, zorder=20)
