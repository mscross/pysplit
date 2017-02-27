"""
=============================================
Basic Scatter Plotting and Advanced MapDesign
=============================================

In ``basic_plotting_example.py``, we learned how to quickly initialize
matplotlib Basemaps with the ``MapDesign`` class and how to plot ``Trajectory``
paths and color according to a single value.  Here, we'll get into more
advanced usage of ``MapDesign`` and learn a couple of ways to scatter plot our
along-trajectory data.

For this example we'll initialize only the February trajectories created in 
``reversetraj_clippedtraj_gen.py``.

"""

import matplotlib.pyplot as plt
import numpy as np

import pysplit

trajgroup = pysplit.make_trajectorygroup(r'C:/trajectories/umn_example/*feb*')

"""
Basemaps and Advanced MapDesign
-------------------------------
We reiterate here that PySPLIT can plot on any matplotlib Basemap, but the
``MapDesign`` class of PySPLIT can expedite setting up Basemaps.

A basic cylindrical map using ``MapDesign`` only requires ``mapcorners`` and
``standard_pm=None``.  However, we can adjust the map projection,
label spacing, arrangement, and sizes; how and which boundaries are drawn, even
switch up the background colors between three different gray level sets,
among other actions.  Basemaps can be drawn on new or existing axes.

Let's make a figure with two Lambert Conformal Conic maps with states and
medium gray backgrounds.  We'll initialize two ``MapDesign``s, since we want
slightly different maps (longitude labels), and make one Basemap per axis.
"""

fig, (ax0, ax1) = plt.subplots(nrows=2, figsize=(8, 10))

mapcorners = [-155, 10, -50, 55]
standard_pm = [-110, 20, 40, 30]

# Note we draw every 15 degrees of longitude, but only label every 30
# Check the pysplit.MapDesign() docs for more options than param_dict shows
param_dict = {'projection':'lcc', 'latlon_labelspacing':(10,30),
              'latlon_spacing':(10,15), 'latlon_fs':16, 'drawstates':True,
              'resolution':'l', 'mapcolor':'medium'}

map_params0 = pysplit.MapDesign(mapcorners, standard_pm, **param_dict)
map_params1 = pysplit.MapDesign(mapcorners, standard_pm, lon_labels=['bottom'], **param_dict)

scattermap0 = map_params0.make_basemap(ax=ax0)
scattermap1 = map_params1.make_basemap(ax=ax1)

"""
Plotting ``Trajectory`` data
----------------------------
This is at present the best way to approximate a color-graded line.
Here, we'll scatter plot along-trajectory pressure data.

Scatter plotting can be performed in a couple ways.  The first is to use the
standard ``basemap.scatter()``, and the second is to use
pysplit.traj_scatter()``.  Both are shown in the example below,
and achieve exactly the same result- the top and bottom maps should be
identical.

There are also different ways of acquiring the ``Trajectory`` coordinates,
which can be mixed and matched with the method of scatter plotting.
"""

# Need to find minimum and maximum pressure first to ensure consistent color
# scaling across all trajectories!
min_pressure = 1000.0
max_pressure = 0.0

for traj in trajgroup:
    if traj.data.Pressure.max() > max_pressure:
        max_pressure = traj.data.Pressure.max()
    if traj.data.Pressure.min() < min_pressure:
        min_pressure = traj.data.Pressure.min()

# To cut down on crowding, let's plot every other trajectory
for traj in trajgroup[::2]:
    mappable = scattermap0.scatter(*traj.path.xy,
    								c=traj.data.Pressure.astype(np.float64).values,
    								cmap=plt.cm.magma, vmin=min_pressure,
    								vmax=max_pressure, zorder=15, latlon=True,
    								edgecolor='none')
    
    mappable = pysplit.traj_scatter(traj.data.Pressure.astype(np.float64).values,
                                    traj.data.geometry.apply(lambda p: p.x).values,
                                    traj.data.geometry.apply(lambda p: p.y).values,
                                    scattermap1, colormap=plt.cm.magma,
                                    vmin=min_pressure, vmax=max_pressure,
                                    suppress_printmsg=True)

"""
Another Look at ``pysplit.traj_scatter()``
-------------------------------------------

A potential advantage of this method is that is can streamline the
application of various normalizations.  Here's a BoundaryNorm dividing our
along-trajectory pressures up into six levels.

Later examples will show other usages of ``pysplit.traj_scatter()``,
colorbar generation, and how to plot ``Trajectory.uptake()`` data.
"""

map_params2 = pysplit.MapDesign(mapcorners, standard_pm, **param_dict)
boundarymap = map_params2.make_basemap()

for traj in trajgroup[::2]:    
    mappable = pysplit.traj_scatter(traj.data.Pressure.astype(np.float64).values,
                                    *traj.path.xy, hymap=boundarymap,
                                    colormap=plt.cm.magma, vmin=min_pressure,
                                    vmax=max_pressure, suppress_printmsg=True,
                                    cnormalize='boundary', levels=6)
