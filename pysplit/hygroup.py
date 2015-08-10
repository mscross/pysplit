class HyGroup(object):
    """
    Superclass for TrajectoryGroups and ClusterGroups

    """

    def __init__(self, trajectories):
        """
        Initialize ``HyPath``.

        Parameters
        ----------
        trajectories : list of ``Trajectory`` instances
            ``Trajectory`` instances that belong to the group

        """

        self.trajectories = trajectories
        self.trajcount = len(trajectories)

    def addgroups(self, other):
        """
        Create new group from the set of ``Trajectory`` instances
        in two groups.

        Parameters
        ----------
        other : ``HyGroup`` or ``HyGroup`` subclass instance

        Returns
        -------
        newgroup : list
            A list of the unique trajectories from the two groups

        """

        set0 = set(self.trajectories)
        set1 = set(other.trajectories)

        newgroup = list(set0 | set1)

        return newgroup
