from __future__ import division, print_function

from .hygroup import HyGroup


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
        Get ``Trajectory`` or ``TrajectoryGroup``.

        Parameters
        ----------
        index : int or slice

        Returns
        -------
        ``Trajectory`` or ``TrajectoryGroup`` depending if indexed
        or sliced.  Won't return a ``Cluster`` because those are
        specially defined.

        """
        newthing = self.trajectories[index]

        if isinstance(newthing, list):
            newthing = TrajectoryGroup(newthing)

        return newthing

    def __add__(self, other):
        """
        Add a ``HyGroup`` to this ``TrajectoryGroup`` instance.

        Parameters
        ----------
        other : ``HyGroup``
            Another ``TrajectoryGroup`` or ``Cluster``.  May or may not
            contain some of the same ``Trajectory`` instances

        Returns
        -------
        A new ``TrajectoryGroup`` containing the union of the sets
        of ``Trajectory`` instances.

        """
        return TrajectoryGroup(HyGroup.__add__(self, other))

    def __sub__(self, other):
        """
        Subtract a ``HyGroup`` from this ``TrajectoryGroup`` instance.

        Parameters
        ----------
        other : ``HyGroup``
            Another ``TrajectoryGroup`` or ``Cluster``

        Returns
        -------
        A new ``TrajectoryGroup`` containing the set difference between
        the sets of ``Trajectory`` instances.

        """
        return TrajectoryGroup(HyGroup.__sub__(self, other))

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
