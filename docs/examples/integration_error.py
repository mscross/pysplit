"""
============================
Estimating Integration Error
============================

Objectives
----------

* Explain purpose of integration error
* Demonstrate integration error calculation procedure with PySPLIT
* Demonstrate another ``TrajectoryGroup``, ``Trajectory`` workflow tool

Intro
-----

Total trajectory error consists of physical and numerical error
(see http://www.arl.noaa.gov/faq_hg11.php).  One part of the numerical error is
integration error, which can be estimated by generating an original/reverse
trajectory pair, in which the reverse trajectory is initialized at the end of
the original trajectory and run the opposite direction.

By calculating the total travel distance of the two trajectories and the
distance between the original trajectory start and reverse trajectory end
points, we can estimate absolute and relative integration error.

First, reverse trajectories must be available.  For information on how to
generate reverse trajectories with ``PySPLIT``, see ``bulk_trajgen_example.py``
and ``reversetraj_clippedtraj_gen.py``.

Setup
-----

Load the original and reverse trajectories.  This example uses the
trajectories generated in ``bulk_trajgen_example.py``.

"""
from __future__ import print_function
import numpy as np
import pysplit

trajgroup = pysplit.make_trajectorygroup(r'C:/trajectories/colgate/*')

for traj in trajgroup:
    traj.load_reversetraj()
    
"""
Calculating integration error
-----------------------------

Values computed when calling ``Trajectory.calculate_integrationerr()``:
    ``Trajectory.integration_error``, the relative error (%)
    ``Trajectory.integration_error_abs``, the absolute error (meters)

"""
for traj in trajgroup:
    traj.calculate_integrationerr()

"""
Usage example
-------------

Once we can have these values, one action we can take is to discard the "bad"
trajectories.  A reasonable way to define "bad" trajectories are those with
integration errors greater than two standard deviations above the mean:

"""
relative_errors = [traj.integration_error for traj in trajgroup]
cutoff = np.mean(relative_errors) + (np.std(relative_errors) * 2)

print('Integration error upper limit: ', cutoff)

"""
With this data, we can cycle through ``trajgroup`` and either identify "good"
trajectories to put in a new ``TrajectoryGroup``
(see ``traj_trajgroup_basics_example.py``), or, as below, identify "bad"
trajectories to remove from ``trajgroup``.  In this example, we make a list of
the identifiers (``Trajectory.trajid``) of "bad" trajectories, then pass the
list to the ``TrajectoryGroup.pop()`` method, which removes the indicated
trajectories from ``trajgroup``.

``TrajectoryGroup.pop()`` accepts a list of ``Trajectory.trajid``, a single 
``Trajectory.trajid``, or an index.  If none of the above are specified it
defaults to the last ``Trajectory`` in ``TrajectoryGroup``.  As with
``list.pop()``, performing ``TrajectoryGroup.pop()``while iterating over
``TrajectoryGroup`` will lead to unexpected behavior. 

``TrajectoryGroup.pop()`` returns a ``Trajectory``, if one
``Trajectory.trajid`` or an index is given, or a ``TrajectoryGroup``,
if given a list of ``Trajectory.trajid``.

"""
bad = []
for traj in trajgroup:
    if traj.integration_error > cutoff:
        bad.append(traj.trajid)

print('Expectation: ', trajgroup.trajcount, 'trajectories -', len(bad),
      'bad trajectories =', trajgroup.trajcount-len(bad), 'trajectories')

trajgroup.pop(trajid=bad)

print('Result: ', trajgroup.trajcount, 'trajectories')

# for traj in trajgroup:
#     print(traj.integration_error)
