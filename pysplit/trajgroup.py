from __future__ import division, print_function
import os
from hygroup import HyGroup


class TrajectoryGroup(HyGroup):
    """
    Class for processing and plotting multiple ``Trajectory`` instances.

    :subclass: of ``HyGroup``.

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

    def pop(self, ind=-1, trajid=None):
        """
        Intercept ``HyGroup.pop()``.

        If a list of ``Trajectory`` instances is returned from
        ``HyGroup.pop()``, return a new ``TrajectoryGroup``.

        Parameters
        ----------
        ind : int
            Default -1.  The positional argument of the ``Trajectory``
            to remove.
        trajid : string
            Default None.  The named argument of the ``Trajectory``
            to remove.  Overrides ``ind`` if not None.

        Returns
        -------
        popped : ``Trajectory`` instance or ``TrajectoryGroup``
            The indicated ``Trajectory`` or a ``TrajectoryGroup`` if multiple
            trajectories were popped out.  May also be a ``None``
            if no matching ``Trajectory`` instances were found.

        """
        popped = HyGroup.pop(self, ind, trajid)

        try:
            popped = TrajectoryGroup(popped)
        except:
            pass

        return popped
