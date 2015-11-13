# PySPLIT

A package for generating [HYSPLIT] (http://ready.arl.noaa.gov/HYSPLIT.php) air parcel trajectories trajectories, performing moisture uptake analyses, expediting HYSPLIT cluster analysis, and for visualizing trajectories, clusters, and along-trajectory meteorological data.  For a basic overview of PySPLIT, see the [SciPy 2015 conference proceedings] (http://conference.scipy.org/proceedings/scipy2015/mellissa_cross_p.html).

## Recent Updates

* PySPLIT now uses the power of GeoPandas rather than pure NumPy
* Faster trajectory file loading/``Trajectory`` object initialization
* The class structure of PySPLIT has been rewritten:
  * ``Trajectory`` and ``Cluster`` objects are now subclasses of ``HyPath`` class, which in turn is a [GeoPandas] (http://geopandas.org/) ``GeoDataFrame`` subclass.
  * ``TrajectoryGroup`` and ``Cluster`` classes are now subclasses of the ``HyGroup`` class.
  * ``HyPath`` and ``HyGroup`` are only used internally, so the API remains essentially the same.
* ``Cluster`` objects are no longer iterable over their member ``Trajectory`` objects (instead iterate over ``Cluster.trajectories``)
* Trajectory generator updates:
  * Improved efficiency
  * Improved API
  * Use *any* weekly or semi-monthly meteorology data (see docs for required filename format), not just gdas1
  * Generate trajectories for every day in each month *OR* for particular slice of days in each month
  
## Using PySPLIT

Coming soon.  Example code snippets are available in the [SciPy 2015 conference proceedings] (http://conference.scipy.org/proceedings/scipy2015/mellissa_cross_p.html).  Much of the work done to PySPLIT since has been under the hood, so these snippets are still a good place to start.  Please cite the proceedings if you use PySPLIT in your work!
