"""
===========================================================
Getting Started with ``TrajectoryGroup`` and ``Trajectory``
===========================================================

The basic PySPLIT workflow consists of cycling through ``TrajectoryGroup``
containers of ``Trajectory`` objects and acting on each ``Trajectory``,
refining the ``TrajectoryGroup`` or creating new ones as necessary.

Initializing the first ``TrajectoryGroup``
------------------------------------------
The first ``TrajectoryGroup`` requires ``Trajectory`` objects to be
initialized from trajectory files.  Here we initialize all of the trajectories
created in ``bulk_trajgen_example.py``.

"""
from __future__ import print_function

import pysplit

trajgroup = pysplit.make_trajectorygroup(r'C:/trajectories/colgate/*')

"""
Workflow
--------
Cycle through the ``TrajectoryGroup`` and act on each trajectory.  Below 
are some sample geometry calculations.

"""
for traj in trajgroup:
    traj.calculate_distance()
    traj.calculate_vector()

"""
Let's create a new ``TrajectoryGroup`` with a list of some of the
``Trajectory`` objects in ``tg``.  For this example, we'll make
a group consisting only of trajectories with rainfall at timepoint 0.

Note- this will only work as written if rainfall was selected as
an output meteorological variable during trajectory generation.  Alternatively,
``Trajectory.set_rainstatus()`` can use humidity variables if those
are available and if the appropriate kwargs are provided.

"""
rainylist = []

for traj in trajgroup:
    traj.set_rainstatus()
    if traj.rainy:
        rainylist.append(traj)
        
rainy_trajgroup = pysplit.TrajectoryGroup(rainylist)

"""
A new ``TrajectoryGroup`` can also be created by addition or subtraction.

Let's subtract our new ``rainy_trajgroup`` from the original
``TrajectoryGroup``, ``trajgroup``.  This yields a new ``TrajectoryGroup``
with only non-rainfall producing trajectories.  The number of 
member ``Trajectory`` objects can be checked using
``TrajectoryGroup.trajcount``

"""
dry_trajgroup = trajgroup - rainy_trajgroup

print(dry_trajgroup.trajcount)
print(rainy_trajgroup.trajcount)
print(trajgroup.trajcount)
