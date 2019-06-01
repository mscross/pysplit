"""
PySPLIT package containing tools for automatically
generating trajectories, performing moisture uptake analyses,
enhancing the HYSPLIT cluster analysis experience,
and visualizing trajectories, trajectory clusters, and
meteorological data along trajectories.

"""

__all__ = ['Trajectory',
           'TrajectoryGroup',
           'print_clusterprocedure',
           'Cluster',
           'ClusterGroup',
           'MapDesign',
           'traj_scatter',
           'traj_path',
           'meteo_contouring',
           'adjust_contourparams',
           'make_cbar',
           'make_cax_cbar',
           'edit_cbar',
           'map_labeller',
           'labelfile_reader',
           'labelfile_generator',
           'make_trajectorygroup',
           'spawn_clusters',
           'hysplit_filelister',
           'load_hysplitfile',
           'load_clusteringresults',
           'generate_bulktraj']

__version__ = '0.3.5'

from .traj import Trajectory

from .trajgroup import TrajectoryGroup

from .clusgroup import print_clusterprocedure, Cluster, ClusterGroup

from .mapdesigner import MapDesign

from .mapmaker import (traj_scatter, meteo_contouring,
                       adjust_contourparams, make_cbar, make_cax_cbar,
                       edit_cbar)

from .maplabeller import map_labeller, labelfile_reader, labelfile_generator

from .hy_processor import make_trajectorygroup, spawn_clusters

from .hyfile_handler import (hysplit_filelister, load_hysplitfile,
                             load_clusteringresults)

from .trajectory_generator import (generate_bulktraj)
