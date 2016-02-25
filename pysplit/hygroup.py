import os


class HyGroup(object):
    """
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

    def make_infile(self, infile_dir):
        """
        Write ``Trajectory`` paths to a file (INFILE).

        Take the member ``Trajectory`` instances and write
        path to INFILE, used by ``HYSPLIT`` to perform cluster analysis.
        Can also use this list of files to reinitialize
        ``TrajectoryGroup``.

        Parameters
        ----------
        infile_dir : string
            The directory in which to create INFILE

        """
        with open(os.path.join(infile_dir, 'INFILE'), 'w') as infile:

            for traj in self:
                try:
                    output = traj.cfullpath
                except:
                    output = traj.fullpath

                output = output.replace('\\', '/')
                infile.writelines(output + '\n')
                infile.flush()
