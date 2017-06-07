from __future__ import division, print_function

import os
import fnmatch
import numpy as np
import pandas as pd


def hysplit_filelister(signature):
    """
    List all HYSPLIT files matching a given signature.

    Parameters
    ----------
    signature : string
        Signature shared by group of HYSPLIT simulation files from a single or
        multiple model runs (if multiple, must contain same output variables).
        This is a Bash-style signature, not a real expression.  The '*' char is
        a wildcard.  Can include an absolute or relative path, or no path.

    Returns
    -------
    matching_files : list of strings
        List of files matching ``signature``

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

    if len(matching_files) == 0:
        raise LookupError("Could not find any files matching the provided "
                          "signature `{0}`, please check your paths and "
                          "try again.".format(signature))

    return matching_files


def load_hysplitfile(filename):
    """
    Load data from each trajectory into a ``NumPy ndarray``.

    Parameters
    ----------
    filename : string
        The name of a trajectory file

    Returns
    -------
    hydata : (M, N) ndarray of floats or list of ndarrays
        Ndarray with M time steps and N variables representing one trajectory.
        If there are multiple trajectories in a file, ``multiple_traj`` will
        be ``True`` and ``hydata`` will be a list of ndarrays, potentially of
        different sizes.
    pathdata : (M, 3) ndarray of floats or list of ndarrays
        The path information in lon, lat, z.  If there are multiple
        trajectories in a file, ``multiple_traj`` will
        be ``True`` and ``pathdata`` will be a list of ndarrays.
    header : list of N strings
        The column headers for ``hydata`` arrays.  Used to parse ``hydata``
        into different trajectory attributes
    datetime : DateTime index of length M
    multiple_traj : Boolean

    """
    # Every header- first part
    header = ['Parcel Number',
              'Timestep']

    with open(filename, 'r') as hyfile:

        contents = hyfile.readlines()

        skip = False
        atdata = False

        # Entire contents because otherwise it misses last line
        for ind, line in enumerate(contents):
            if skip:
                skip = False
                continue

            if atdata:
                data = [float(x) for x in line.split()]
                if multiline:
                    data.extend([float(x) for x in contents[ind + 1].split()])
                    skip = True

                # Don't need date/timestep info in each file line
                # Got the date, time direction from the OMEGA line
                del data[1:8]
                hydata[arr_ind, :] = data[:2] + data[5:]
                pathdata[arr_ind, :] = data[2:5]
                arr_ind += 1
                continue

            # OMEGA happens first
            if 'OMEGA' in line:
                num_parcels = int(line.split()[0])

                multiple_traj = False
                if num_parcels > 1:
                    multiple_traj = True

                # Number of data rows = length of contents minus the number of
                # lines before OMEGA, OMEGA, andbetween OMEGA and first time pt
                flen = len(contents) - (2 + num_parcels) - ind

                # Determine timepoint freq, get all the initialization dates
                date0 = []
                freq = 'H'

                if 'BACKWARD' in line:
                    freq = '-1H'

                for i in range(1, num_parcels + 1):
                    date0.append(contents[ind + i].split()[:4])

                continue

            # PRESSURE happens second
            if 'PRESSURE' in line:
                new_header = line.split()[1:]
                columns = 12 + len(new_header)
                header.extend(new_header)

                multiline = False
                # Pre-Jan 2017 files may have timepoints running onto second
                # lines, which are always short
                if len(contents[ind + 1]) > len(contents[ind + 2]):
                    multiline = True

                    # Data file is only half as many lines as it looks
                    flen /= 2

                # Initialize empty data array
                # 10 dropped columns include date, time step, etc
                hydata = np.empty((int(flen), columns - 10))
                pathdata = np.empty((int(flen), 3))
                atdata = True
                arr_ind = 0
                continue

    # Catch the vast majority of non-HYSPLIT files if passed
    # Works because the above conditionals fall through; vars never defined
    if 'multiline' not in locals() or 'date0' not in locals():
        raise IOError("The file, `{0}`, does not appear to be "
                      "a valid HYSPLIT file. Please double check "
                      "your paths.".format(filename))

    datestrings = []
    for d in date0:
        datestrings.append("20" +
            "{0:02}{1:02}{2:02}{3:02}".format(*[int(x) for x in d]) +
            '0000')

    # Get pathdata in x, y, z from lats (y), lons (x), z
    pathdata = pathdata[:, np.array([1, 0, 2])]

    # Split hydata into individual trajectories (in case there are multiple)
    if multiple_traj:
        hydata, pathdata, datetime = trajsplit(hydata, pathdata, datestrings,
                                               freq)
    else:
        datetime = pd.date_range(datestrings[0], freq=freq,
                                 periods=pathdata.shape[0])

    return hydata, pathdata, header, datetime, multiple_traj


def trajsplit(hydata, pathdata, datestrings, freq):
    """
    Split an array of hysplit data into list of unique trajectory arrays.

    Parameters
    ----------
    hydata : (L, N) ndarray of floats
        Array with L rows and N variables, introspected from a hysplit
        data file.
    pathdata : (L, 3) ndarray of floats
        Array with L rows and x, y z (lons, lats, altitude) columns
    datestrings : list of strings
        String representations of trajectory initialization date and time
    freq : string
        Either 'H' or '-1H', indicating a forwards or backwards trajectory

    Returns
    -------
    split_hydata : list of (?, N) ndarrays of floats
        ``hydata`` split into individual trajectories
    split_pathdata : list of (?, 3) ndarrays of floats
        ``pathdata split into individual trajectories
    datetime : list of pandas date ranges
        List of ranges of DateTime objects for each array within
        ``split_hydata`` and ``split_pathdata``
    """
    # Find number of unique trajectories within `hydata`
    unique_traj = np.unique(hydata[:, 0])

    # Sort the array row-wise by the first column
    # Timepoints from same traj now grouped together
    sorted_indices = np.argsort(hydata[:, 0], kind='mergesort')
    sorted_hydata = hydata[sorted_indices, :]
    sorted_pathdata = pathdata[sorted_indices, :]

    # Find first occurrence of each traj, except for the first
    # which is obviously 0
    first_occurrence = [np.nonzero(sorted_indices == u)[0][0]
                        for u in unique_traj[1:]]

    # Split `hydata` and `pathdata` into list of arrays, one
    # array per traj.  May or may not be equal sizes
    split_hydata = np.split(sorted_hydata, first_occurrence)
    split_pathdata = np.split(sorted_pathdata, first_occurrence)

    datetime = []
    for p, d in zip(split_pathdata, datestrings):
        datetime.append(pd.date_range(d, freq=freq, periods=p.shape[0]))

    return split_hydata, split_pathdata, datetime


def load_clusteringresults(clusterfilename):
    """
    Load a 'CLUSLIST_#' file into an array.

    The 'CLUSLIST_#' file contains the information on how trajectories
    within a trajectory group are distributed among clusters.  Called by
    ``hy_processor.spawn_clusters()``

    Parameters
    ----------
    clusterfilename : string
        Path to CLUSLIST_# file

    Returns
    -------
    traj_inds : list of lists of ints
        Lists of lists of trajectory indices.
    totalclusters : int
        Number of unique clusters.  Corresponds to number of arrays in
        ``clusterarray_list``, number of sublists in ``traj_inds``.

    """
    orig_dir = os.getcwd()

    try:
        head, tail = os.path.split(clusterfilename)

        os.chdir(head)

        with open(tail, 'r') as clusterfile:
            contents = clusterfile.readlines()
            clusterinfo = np.empty((len(contents)))
            traj_inds = np.empty((len(contents)))

            for ind, line in enumerate(contents):

                data = [int(x) for x in line.split()[:-1]]
                clusterinfo[ind] = data[0]
                traj_inds[ind] = data[-1]

        uniqueclusters = np.unique(clusterinfo)
        totalclusters = int(np.max(clusterinfo))

        # Fix off by one in trajectory indicies
        traj_inds = traj_inds - 1

        # Get the indices of the first occurrence of each unique cluster number
        first_occurrence = [np.nonzero(clusterinfo == u)[0][0]
                            for u in uniqueclusters[1:]]

        # Split into separate arrays at the first occurrences
        cluster_trajlist = np.split(traj_inds, first_occurrence)

    finally:
        os.chdir(orig_dir)

    return cluster_trajlist, totalclusters
