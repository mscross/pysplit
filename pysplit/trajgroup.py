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
        Remove Trajectory object(s) from self.

        Shortcut to self.trajectories.pop() that updates the
        self.trajcount and the list of trajids.

        Parameters
        ----------
        ind : int
            The positional argument of the ``Trajectory``
            to remove.
        trajid : string or list of strings
            The identifier(s) of the ``Trajectory`` object(s)
            to remove from ``self``.  Overrides ``ind`` if not None.

        Returns
        -------
        popped : ``Trajectory`` or ``TrajectoryGroup``
            A``Trajectory`` or ``TrajectoryGroup`` consisting of the
            trajectory or trajectories indicated by ``ind`` or ``trajid``.

        """
        if trajid is not None:
            try:
                to_pop = [self.trajids.index(trajid)]
            except ValueError:
                to_pop = [self.trajids.index(t) for t in trajid
                          if t in self.trajids]
                if len(to_pop) == 0:
                    raise ValueError('TrajIDs not in list of self.trajids')

            to_pop.sort()
            popped = []
            for p in to_pop[::-1]:
                popped.append(self.trajectories.pop(p))
                self.trajids.pop(p)
            self.trajcount = len(self.trajectories)

            if len(popped) == 1:
                popped = popped[0]
            else:
                popped = TrajectoryGroup(popped)
        else:
            popped = self.trajectories.pop(ind)
            self.trajids.pop(ind)
            self.trajcount = len(self.trajectories)

        return popped

    def append(self, traj):
        """
        Add a ``Trajectory`` to the ``self``.

        Parameters
        ----------
        traj : ``Trajectory`` instance
            The ``Trajectory`` to add to the end of ``self``.

        """
        if hasattr(traj, 'trajid'):
            self.trajectories.append(traj)
            self.trajids.append(traj.trajid)
            self.trajcount = len(self.trajectories)
