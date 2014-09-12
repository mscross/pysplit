import numpy as np
import math
import os
import fnmatch


def hysplit_filelister(signature):
    """
    List all HYSPLIT files matching a given signature.

    Parameters
    ----------
    signature : string
        Signature shared by group of HYSPLIT simulation files from a single or
        multiple model runs (if multiple, must contain same output variables).
        This is a Bash-style signature, not a real expression.  The `*` char is
        a wildcard.  Can include an absolute or relative path, or no path.

    Returns
    -------
    matching_files : list of strings
        List of files matching `signature`

    Notes
    -----
    Any Bash-style signature is supported.
        The file search is non-recursive.

    """

    # Initialize
    orig_dir = os.getcwd()
    matching_files = []

    try:
        head, tail = os.path.split(signature)

        os.chdir(head)

        # os.walk obtains list in top-down manner (files in order)
        _, _, files = next(os.walk('.'))

        for each_file in files:
            if fnmatch.fnmatch(each_file, tail):
                matching_files.append(each_file)

    finally:
        os.chdir(orig_dir)


    return matching_files



def load_hysplitfile(filename):
    """
    Load data from each trajectory (a single hysplit file) into a NumPy ndarray.

    Parameters
    ----------
    filename : string
        The name of a trajectory file

    Returns
    -------
    hydata : list of (M, N) ndarrays of floats
        List of arrays with M time steps and N variables.  Each array
        represents one trajectory; typically but not always there is one
        trajectory per file
    header : list of N strings
        The column headers for `hydata` arrays.  Used to parse hydata into
        different trajectory attributes
    filelist : list of identical strings
        len(filelist)= len(hydata)
        The filename.  Eventually becomes attribute in trajectory object

    """

    hyfile = open(filename, 'r')

    # Get through meteofile info at the head of the file
    while True:
        line = hyfile.readline()
        if 'RAINFALL' in line:
            break

    # Begin constructing header
    header = ['Parcel Number',
              'Year',
              'Month',
              'Date',
              'Hour (UTC)',
              '24 mod 6 (hr)',
              'Time step (hr)',
              'Latitude',
              'Longitude',
              'Altitude (magl)']

    # Acquire rest of header
    new_header = line.split()[1:]
    columns = 10 + 2 + len(new_header)
    header.extend(new_header)

    # Check if data is multiline
    # Multiline data not supported for some applications
    multiline = False
    if columns > 20:
        multiline = True

    # Initialize empty data array
    hydata = np.empty((0, columns - 2))

    while True:
        line = hyfile.readline()
        if line =='':
            break

        time_step = line.split()

        # Check if data is multiline.  Multi-line data not supported for
        # all applications
        if multiline:
            time_step.extend(hyfile.readline().split())

        # Put into array with unnecessary data sliced out, concatenate
        time_step = np.asarray(time_step).astype(np.float64)
        time_step = np.r_[time_step[0], time_step[2:6], time_step[7:]]
        time_step = np.atleast_2d(time_step)

        if time_step.shape[1] != hydata.shape[1]:
            raise ValueError("Introspected number of columns incorrect!")

        hydata = np.concatenate((hydata, time_step), axis=0)

    # Split hydata into individual trajectories (in case there are multiple)
    hydata = trajsplit(hydata)
    filelist = [filename] * len(hydata)

    return hydata, header, filelist



def trajsplit(hydata):
    """
    Splits an array of hysplit data into list of unique trajectory arrays

    Parameters
    ----------
    hydata : (L, N) ndarray of floats
        Array with L rows and N variables, introspected from a hysplit
        data file.

    Keyword Arguments
    -----------------

    Returns
    -------
    split_hydata : list of (M, N) ndarrays of floats
        `hydata` split into individual trajectories
    OR
    [hydata] : `hydata` in list

    """

    # Find number of unique trajectories within `hydata`
    unique_traj = np.unique(hydata[:,0])
    if unique_traj.size > 1:

        # Sort the array row-wise by the first column
        traj_id = hydata[:,0]
        sorted_indices = np.argsort(traj_id, kind='mergesort')
        sorted_hydata = hydata[sorted_indices, :]

        # Split `hydata` into a list of arrays of three equal sizes
        split_hydata = np.split(sorted_hydata, unique_traj.size)

        return split_hydata

    else:

        return [hydata]



def load_clusterfile(clusterfilename):
    """
    Load a CLUSLIST_# file into an array

    The CLUSLIST_# file contains the information on how trajectories
        within a trajectory group are distributed among clusters

    Parameters
    ----------
    clusterfilename : string
        Path to CLUSLIST_# file

    Keyword Arguments
    -----------------

    Returns
    -------
    traj_inds : list of lists of ints
        Lists of lists of trajectory numbers belonging to each cluster
    totalclusters : int
        Number of unique clusters.  Corresponds to number of arrays in
        `clusterarray_list`, number of sublists in `traj_inds`.

    Notes
    -----
    `load_clusterfile` returns a list of trajectory indices that may be
        used to pull trajectories out of a TrajectoryGroup's list to create
        Cluster objects and a ClusterGroup.
    All trajectories in a TrajectoryGroup are used in the clustering
        process, so if a subset of the original TrajectoryGroup were desired,
        a new TrajectoryGroup should be created with only that subset, and
        the cluster process should be performed on the new TrajectoryGroup.

    """

    orig_dir = os.getcwd()
    clusterfiledata = []

    try:
        head, tail = os.path.split(clusterfilename)

        os.chdir(head)

        clusterfile = open(tail, 'r')

        # Read in clustering information
        while True:
            line = clusterfile.readline()
            if line == '':
                break

            dataline = line.split()

            # Collect cluster number and trajectory number
            dataline = [dataline[0]] + [dataline[-2]]
            clusterfiledata.append(dataline)

        # Stack lines into an array
        filearray = np.vstack(clusterfiledata)

        clusternums = filearray[:,0].astype(int)

        # Find unique cluster numbers and their total number
        uniqueclusters = np.unique(clusternums)
        totalclusters = np.max(clusternums)

        # Get the indices of the first occurrence of each unique cluster number
        first_occurrence = []
        for u in uniqueclusters:
             first_occurrence.append(np.nonzero(clusternums == u)[0][0])

        # Split into separate arrays at the first occurrences
        # First occurrence of first cluster is not considered,
        # since it is zero and will result in an empty array
        clusterarray_list = np.vsplit(filearray, first_occurrence[1:])

        # For each cluster, create a list of indices corresponding to the
        # constituent trajectories
        traj_inds = []
        for cluster in clusterarray_list:
            inds = [int(i)-1 for i in cluster[:,1]]
            traj_inds.append(inds)

    finally:
        os.chdir(orig_dir)

    return traj_inds, totalclusters