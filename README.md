# PySPLIT

A package for generating [HYSPLIT](http://ready.arl.noaa.gov/HYSPLIT.php) air parcel trajectories trajectories, performing moisture uptake analyses, expediting HYSPLIT cluster analysis, and for visualizing trajectories, clusters, and along-trajectory meteorological data.

For an overview and brief history of PySPLIT, a new, updated technical paper- Introduction to PySPLIT: A Python Toolkit for NOAA ARL’s HYSPLIT Model- can be found in the Sept/Oct 2018 (vol 20, issue 5, p. 47-62 ) issue of *Computing in Science and Engineering*!  This supercedes the [SciPy 2015 conference proceedings](http://conference.scipy.org/proceedings/scipy2015/mellissa_cross_p.html).

### If you are running version 0.3.3 or older and are performing moisture uptake analyses, please update PySPLIT and rerun your moisture uptake analyses immediately.  Geographic points were previously assigned to ``Trajectory.uptake`` backwards.  This has been corrected.

## Coming Soon
* HYSPLIT clustering fully in PySPLIT
* Increased trajectory generation functionality:
  * New modes
  * More control over trajectory initialization conditons
  * Improved meteorology discovery and better support for sub-weekly files
* Support for matrix and ensemble trajectories
* Extended library of examples
* Various quality of life/convenience updates and more!


## Past Updates

* Support for Cartopy (basemap use to be deprecated in future update)
* Support for Python 3.6 and 3.7
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
  * Use *any* weekly or semi-monthly meteorology data (see docs for required filename format), not just gdas1, and not just from the 21st century!
  * Generate trajectories for every day in each month *OR* for particular slice of days in each month
  * Generate reverse trajectories at time of bulk trajectory generation OR during your analysis workflow!
* Choose the starting point for your moisture uptake analyses
* Removal of certain assumptions about trajectory file structure for HYSPLIT January 2017 (854) compatibility.
* Removal of certain assumptions about trajectory century and timepoint interval from loading process
* Check out the growing library of examples!
## Installing PySPLIT

PySPLIT is compatible with Python 2.7, 3.6, and 3.7.  It depends on:
* NumPy >= 1.6
* matplotlib >= 1.2
* Basemap >= 1.0
* GeoPandas >= 0.1
* Cartopy >= 0.15

and is available on PyPi.  You can install the latest stable release by running:

```
$ pip install pysplit
```

To install from source or create a development installation, clone and fork PySPLIT then install by running:

```
$ python setup.py install
```

or develop locally by running:

```
$ python setup.py develop
```

### Installing in a conda virtual environment:

Installation difficulties with PySPLIT are typically related to GeoPandas dependencies.  An easy work-around is installing PySPLIT in a new conda virtual environment.  This is the recommended installation method.  First, add the conda-forge channel:
```
$ conda config --add channels conda-forge
```

Next, create the conda environment.  For a Python 3.6 environment named `pysplitenv`, run:
```
$ conda create --name pysplitenv python=3.6 numpy matplotlib pandas basemap six fiona shapely geopandas cartopy
```

Similarly, for a Python 3.7 environment named `pysplitenv`, run:
```
$ conda create --name pysplitenv python=3.7 numpy matplotlib pandas basemap six fiona shapely geopandas cartopy
```

Or, to create a Python 2.7 environment named `pysplitenv`, run:
```
$ conda create --name pysplitenv python=2.7 numpy matplotlib pandas basemap six fiona=1.5.1 shapely geopandas cartopy
```

Activate `pysplitenv` by running the following on Windows:
```
$ activate pysplitenv
```
If you are on Linux or OSX, instead run:
```
$ source activate pysplitenv
```

Within your virtual environment, install PySPLIT as above.

## Using PySPLIT

Examples can be found in docs/examples.  PySPLIT is currently tested on Windows 7 using HYSPLIT revision 927 (Feb. 2018) and the preferred PySPLIT installation methods listed above.   Many thanks are due to the NOAA Air Research Laboratory for providing the HYSPLIT model.
