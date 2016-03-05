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
        self.trajids = [traj.trajid for traj in self.trajectories]

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

    def pop(self, ind, trajid):
        """
        Remove Trajectory object from self.

        Shortcut to self.trajectories.pop() that updates the
        self.trajcount and the list of trajids.

        Parameters
        ----------
        ind : int
            The positional argument of the ``Trajectory``
            to remove.
        trajid : string
            The named argument of the ``Trajectory``
            to remove.  Overrides ``ind`` if not None.

        Returns
        -------
        popped : ``Trajectory`` instance or list of
            The indicated ``Trajectory`` or a list if multiple
            trajectories were popped out.  May also be a ``None``
            if no matching ``Trajectory`` instances were found.

        """
        if trajid is not None:
            popped = []
            to_pop = [self.trajids.index(s) for s in
                      self.trajids if trajid in self.trajids]
            for p in to_pop[::-1]:
                popped.append(self.trajectories.pop(p))
                self.trajids.pop(p)
            self.trajcount = len(self.trajectories)
            if len(popped) == 0:
                popped = None
                print('No trajectories found matching ', trajid)
            elif len(popped) == 1:
                popped = popped[0]
        else:
            popped = self.trajectories.pop(ind)
            self.trajids.pop(ind)
            self.trajcount = len(self.trajectories)

        return popped

    def make_infile(self, infile_dir, use_clippedpath=True):
        """
        Take ``Trajectory`` instances in ``HyGroup`` and write
        path to INFILE, used by ``HYSPLIT`` to perform cluster analysis.

        If a specific subset of ``Trajectory`` instances is needed,
        create a new ``HyGroup`` containing only qualifying
        ``Trajectory`` instances.

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
                if use_clippedpath:
                    try:
                        output = traj.cfullpath
                    except:
                        output = traj.fullpath
                else:
                    output = traj.fullpath

                output = output.replace('\\', '/')
                infile.writelines(output + '\n')
                infile.flush()