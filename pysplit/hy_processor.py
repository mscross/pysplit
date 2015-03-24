from __future__ import division
import os
import numpy as np
import hyfile_handler as hh
from traj import Trajectory
from trajgroup import TrajectoryGroup
from clusgroup import Cluster, ClusterGroup


def make_trajectorygroup(signature):
    """
    Initialize trajectory objects from HYSPLIT back trajectory data files.

    Parameters
    ----------
    signature : string
        Signature shared by a group of HYSPLIT simulation files from one or
        multiple model runs (if multiple, must contain same output variables).
        This is a Bash-style signature, not a real expression.  The `*` char is
        a wildcard.  Can include an absolute or relative path, or no path.

    Returns
    -------
    trajectories : trajectory group
        Trajectory group object containing trajectory objects created from all
        simulation files matching `signature`.

    """

    # Get list of hysplit files matching signature
    hyfiles = hh.hysplit_filelister(signature)

    orig_dir = os.getcwd()

    traj_list = []
    filename_list = []
    trajectories = []

    try:
        head, _ = os.path.split(signature)

        os.chdir(head)

        # Sort list of hysplit files by the datestring at the end
        # Will also sort in ascending altitude within each same datestring
        hyfiles.sort(key=lambda x: x[-8:])

        # Load in the hysplit file data
        # Get lists of the datestrings and filenames of the hysplit files
        for hyfile in hyfiles:

            hydata, header, filelist = hh.load_hysplitfile(hyfile)

            traj_list.extend(hydata)
            filename_list.extend(filelist)

        # Initialize trajectory objects
        for data, filename in zip(traj_list, filename_list):

            fullpath = os.path.join(head, filename)

            trajectory = Trajectory(data, header, fullpath)

            trajectories.append(trajectory)

        # initialize trajectory group
        trajectories = TrajectoryGroup(trajectories)
    finally:
        os.chdir(orig_dir)

    return trajectories, filename_list


def spawn_clusters(traj_group, cfile, endpoint_dir):
    """
    Acquires the distribution of trajectories from `cfile` and
        creates new cluster objects and a cluster group object based
        on that information

    Parameters
    ----------
    cfile : string
        The filename of the CLUSLIST_# file that indicates
        trajectory distribution among clusters

    Returns
    -------
    clustergroup : ClusterGroup object
        The group of clusters derived from original TrajectoryGroup.
        A ClusterGroup consists of a list of Cluster objects, which are
        specialized TrajectoryGroup objects.

    """

    traj_inds, totalclusters = hh.load_clusterfile(cfile)

    all_clusters = []

    for i in range(0, totalclusters):
        # Get cluster number and pick out member trajectories
        cluster_num = i + 1
        trajlist = [traj_group.trajectories[j] for j in traj_inds[i]]

        # Get the cluster path
        endpt_fname = ('C' + str(cluster_num) + '_' +
                       str(totalclusters) + 'mean.tdump')
        endpt_file = os.path.join(endpoint_dir, endpt_fname)
        data, header, _ = hh.load_hysplitfile(endpt_file)

        latitude = data[0][:, header.index('Latitude')]
        longitude = data[0][:, header.index('Longitude')]

        # Make sure longitudes are -180 to 180
        longitude = np.where(longitude > 180.0, longitude - 360.0, longitude)

        # Make cluster
        clusterobj = Cluster(trajlist, cluster_num, latitude, longitude)

        all_clusters.append(clusterobj)

    clustergroup = ClusterGroup(all_clusters)

    return clustergroup
