from __future__ import division, print_function
import os
from hygroup import HyGroup


class TrajectoryGroup(HyGroup):
    """
    Class for processing and plotting multiple ``Trajectory`` instances.

    """

    def __init__(self, trajectories):
        """
        Initialize ``TrajectoryGroup`` object.

        Parameters
        ----------
        trajectories : list of ``Trajectory`` instances
            ``Trajectory`` instances that belong in the group.

        """
        HyGroup.__init__(self, trajectories)

    def __getitem__(self, index):
        """
        Index or slice ``self.trajectories`` to get a ``Trajectory`` or
        ``TrajectoryGroup``, respectively

        """

        newthing = self.trajectories[index]

        # TrajectoryGroup requires a list of Trajectory instances,
        # but won't fail if given a single Trajectory
        if isinstance(newthing, list):
            newthing = TrajectoryGroup(newthing)

        return newthing

    def addgroups(self, other):
        """
        Create new ``TrajectoryGroup`` from two ``TrajectoryGroup`` instances.
        Checks for duplicate ``Trajectory`` instances.

        Parameters
        ----------
        other : ``TrajectoryGroup`` or ``Cluster``
            A different ``TrajectoryGroup`` or ``Cluster`` that may or may not
            contain some of the same ``Trajectory`` instances

        Returns
        -------
        new_self : ``TrajectoryGroup``
            A new ``TrajectoryGroup`` from the combination of
            ``self`` and ``other``

        """

        new_tg = TrajectoryGroup(HyGroup.addgroups(self, other))

        return new_tg

    def make_infile(self, infile_dir):
        """
        Take ``Trajectory`` instances in ``TrajectoryGroup`` and write
        path to infile

        If a specific subset of ``Trajectory`` instances is needed,
        create a new ``TrajectoryGroup`` containing only qualifying
        ``Trajectory`` instances

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
