import numpy as np
import os
import hyfile_handler as hh
from .traj import Trajectory, TrajectoryGroup


def hysplit_file_processor(signature):
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
    datestrings = []

    try:
        head, _ = os.path.split(signature)

        os.chdir(head)

        # Sort list of hysplit files by the datestring at the end
        # Will also sort in ascending altitude within each same datestring
        hyfiles.sort(key = lambda x: x[-8:])

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