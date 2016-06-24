# PySPLIT

A package for generating [HYSPLIT] (http://ready.arl.noaa.gov/HYSPLIT.php) air parcel trajectories trajectories, performing moisture uptake analyses, expediting HYSPLIT cluster analysis, and for visualizing trajectories, clusters, and along-trajectory meteorological data.  For a basic overview of PySPLIT, see the [SciPy 2015 conference proceedings] (http://conference.scipy.org/proceedings/scipy2015/mellissa_cross_p.html).

## Recent Updates

* PySPLIT now uses the power of GeoPandas rather than pure NumPy
* Faster trajectory file loading/``Trajectory`` object initialization
* Need help clustering?  ``pysplit.print_clusteringprocedure()``.
* The class structure of PySPLIT has been rewritten:
  * ``Trajectory`` and ``Cluster`` objects are now subclasses of ``HyPath`` class.
  * Along-trajectory data for ``HyPath`` classes lives in the ``data`` attribute, a [GeoPandas] (http://geopandas.org/) ``GeoDataFrame``.
  * ``TrajectoryGroup`` and ``Cluster`` classes are now subclasses of the ``HyGroup`` class.  They are both iterable; they can also be added together or subtracted.
  * ``HyPath`` and ``HyGroup`` are only used internally, so the API remains essentially the same.
* Trajectory generator updates:
  * Improved efficiency
  * Improved API
  * Use *any* weekly or semi-monthly meteorology data (see docs for required filename format), not just gdas1
  * Generate trajectories for every day in each month *OR* for particular slice of days in each month
  
## Using PySPLIT

Coming soon.  Example code snippets are available in the [SciPy 2015 conference proceedings] (http://conference.scipy.org/proceedings/scipy2015/mellissa_cross_p.html).  Much of the work done to PySPLIT since has been under the hood, so these snippets are still a good place to start.  Please cite the proceedings if you use PySPLIT in your work!
