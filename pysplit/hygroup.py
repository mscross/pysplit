from __future__ import division, print_function

import os


class HyGroup(object):
    """
    Class for initializing ``HyPath`` container object.

    :superclass: for ``TrajectoryGroup`` and ``ClusterGroup``.

    """

    def __init__(self, trajectories):
        """
        Initialize ``HyGroup``.

        Parameters
        ----------
        trajectories : list of ``Trajectory`` instances
            ``Trajectory`` instances that belong to the group

        """
        self.trajectories = trajectories
        self.trajcount = len(trajectories)
        self.trajids = [traj.trajid for traj in self.trajectories]

    def __add__(self, other):
        """
        Add ``HyGroup`` instances.

        Create ``HyGroup`` from the union of two sets of ``Trajectory``
        instances.

        Parameters
        ----------
        other : ``HyGroup`` subclass instance

        Returns
        -------
        newgroup : list
            A list of the unique trajectories from the two groups.
            Used to make new ``TrajectoryGroup`` instance.

        """
        set0 = set(self.trajectories)
        set1 = set(other.trajectories)

        newgroup = list(set0 | set1)

        return newgroup

    def __sub__(self, other):
        """
        Subtract ``HyGroup`` instances.

        Create new ``HyGroup`` from the set difference of two
        sets of ``Trajectory`` instances.

        Parameters
        ----------
        other : ``HyGroup`` subclass instance

        Returns
        -------
        newgroup : list
            A list of the set difference of the trajectories.
            Has has all the elements of ``self`` with the
            trajectories of ``other`` removed.  Used to
            make new ``TrajectoryGroup`` instance.

        """
        set0 = set(self.trajectories)
        set1 = set(other.trajectories)

        newgroup = list(set0 - set1)

        return newgroup

    def make_infile(self, infile_dir, use_clippedpath=True):
        """
        Write member ``HyPath`` file paths to INFILE.

        INFILE is used by ``HYSPLIT`` to perform cluster analysis.
        If a specific subset of ``HyPath`` instances is needed,
        create a new ``HyGroup`` containing only qualifying
        ``HyPath`` instances.

        Parameters
        ----------
        infile_dir : string
            The directory in which to create INFILE

        use_clippedpath : Boolean
            Default True. Write out path of clipped trajectory
            rather than original trajectory.

        """
        with open(os.path.join(infile_dir, 'INFILE'), 'w') as infile:

            for traj in self:
                skip_output = False
                if traj.multitraj:
                    try:
                        output = traj.cfullpath
                    except AttributeError:
                        print(traj.trajid, " missing clusterable file")
                        skip_output = True
                else:
                    if use_clippedpath:
                        try:
                            output = traj.cfullpath
                        except AttributeError:
                            output = traj.fullpath
                    else:
                        output = traj.fullpath

                if not skip_output:
                    output = output.replace('\\', '/')
                    infile.writelines(output + '\n')
                    infile.flush()
