from __future__ import division, print_function

import os
import numpy as np

from .hyfile_handler import (load_hysplitfile, load_clusteringresults,
                             hysplit_filelister)
from .traj import Trajectory
from .trajgroup import TrajectoryGroup
from .clusgroup import Cluster, ClusterGroup


def make_trajectorygroup(signature):
    """
    Initialize ``Trajectory`` instances from HYSPLIT trajectory data files.

    Parameters
    ----------
    signature : string or iterable of strings
        Signature shared by a group of HYSPLIT simulation files from one or
        multiple model runs (if multiple, must contain same output variables).
        This is a Bash-style signature, not a real expression.  The `*` char is
        a wildcard. Can include an absolute or relative path, or no path to
        target the current directory. If ``signature`` is not a string, it is
        assumed to be an iterable containing path(s) to specifc HYSPLIT files.

    Returns
    -------
    trajectories : ``TrajectoryGroup``
        ``TrajectoryGroup`` object containing ``Trajectory`` instances created
        from all simulation files matching ``signature``.

    """
    # Get list of hysplit files matching signature
    if not isinstance(signature, str):
        hyfiles = [os.path.split(hyfile)[-1] for hyfile in signature]
        folder, _ = os.path.split(signature[0])
    else:
        hyfiles = hysplit_filelister(signature)
        folder, _ = os.path.split(signature)

    orig_dir = os.getcwd()
    trajectories = []

    # Wrap in try .. finally to ensure current working directory restored
    try:
        os.chdir(folder)

        # Sort list of hysplit files by the datestring at the end
        # Will also sort in ascending altitude within each same datestring
        hyfiles.sort(key=lambda x: x[-8:])

        clipdir = os.path.join(folder, 'clippedtraj')
        if not os.path.isdir(clipdir):
            clipdir = None

        # Load in the hysplit file data
        # Get lists of the datestrings and filenames of the hysplit files
        for hyfile in hyfiles:

            data, path, head, datetime, multitraj = load_hysplitfile(hyfile)

            if multitraj:
                # Initialize trajectory objects
                for d, p, dt in zip(data, path, datetime):

                    # Get rid of parcel number in d
                    # Get rid of parcel #, lat, lon, altitude in head
                    trajectories.append(Trajectory(d, p, dt, head,
                                                   folder, hyfile, clipdir,
                                                   multitraj))

            else:
                trajectories.append(Trajectory(data, path, datetime, head,
                                               folder, hyfile, clipdir,
                                               multitraj))

        # initialize trajectory group
        trajectories = TrajectoryGroup(trajectories)
    finally:
        # Restore current working directory no matter what
        os.chdir(orig_dir)

    return trajectories


def spawn_clusters(trajgroup, assignment_file, cluspath_dir):
    """
    Make ``ClusterGroup`` from ``trajgroup``.

    Acquires the distribution of ``Trajectories`` in ``trajgroup``
    from ``assignment_file`` and creates new ``Cluster`` and ``ClusterGroup``
    instances based on that information

    Parameters
    ----------
    assignment_file : string
        The name of the 'CLUSLIST_#' file that indicates the assignment of each
        ``Trajectory``  to a ``Cluster``

    Returns
    -------
    clustergroup : `ClusterGroup` instance
        A group of `Clusters` derived from original ``TrajectoryGroup``
        (``trajgroup``).  A ``ClusterGroup`` consists of a list of ``Cluster``
        objects, which are specialized ``TrajectoryGroup`` objects.

    """
    traj_inds, totalclusters = load_clusteringresults(assignment_file)

    all_clusters = []

    for i in range(0, totalclusters):
        # Get cluster number and pick out member trajectories
        cluster_num = i + 1
        trajlist = [trajgroup.trajectories[int(j)] for j in traj_inds[i]]

        # Get the cluster path
        cluspath_fname = ('C' + str(cluster_num) + '_' +
                          str(totalclusters) + 'mean.tdump')
        cluspath_file = os.path.join(cluspath_dir, cluspath_fname)
        data, pathdata, header, datetime, _ = load_hysplitfile(cluspath_file)

        # Make sure longitudes are -180 to 180
        pathdata[:, 1] = np.where(pathdata[:, 1] > 180.0,
                                  pathdata[:, 1] - 360.0,
                                  pathdata[:, 1])

        # Make cluster
        clusterobj = Cluster(data, pathdata, datetime, header, trajlist,
                             cluster_num)

        all_clusters.append(clusterobj)

    clustergroup = ClusterGroup(all_clusters)

    return clustergroup
