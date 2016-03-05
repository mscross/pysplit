from __future__ import division, print_function
from trajgroup import TrajectoryGroup
from hypath import HyPath
from hygroup import HyGroup

"""
Procedure
---------
In ``PySPLIT``
1.  Create ``TrajectoryGroup`` with desired set of trajectories
2.  ``TrajectoryGroup.make_infile()``

In ``HYSPLIT``
3.  Trajectory --> Special Runs --> Clustering --> Standard
4.  Adjust clustering parameters, endpoints (trajectory) folder, and
    working folder (where output will be stored)
5.  ``Run cluster analysis`` and determine  and set appropriate
    number of clusters
6.  Assign trajectories to clusters (``Run``)
7.  ``Display Means`` and ``Display Clusters``, ``Quit``

In ``PySPLIT``
8.  ``spawn_clusters(TrajectoryGroup, distribution_file, endpoint_dir)``

This creates a ``ClusterGroup`` populated by ``Cluster``s.

"""


class Cluster(HyPath, HyGroup):
    """
    A special :subclass: of both ``HyGroup`` and ``HyPath``.

    Clusters contain both trajectories and mean path information.  The mean
    path and the trajectory composition is determined by ``HySPLIT``.

    Clusters are not iterable over trajectories in order to avoid
    conflicts with the path data.

    """

    def __init__(self, clusterdata, pathdata, datetime, clusterheader,
                 trajectories, cluster_number):
        """
        Initialize ``Cluster`` object.

        Parameters
        ----------
        trajectories : list of ``Trajectory`` objects
            Trajectories that belong in the cluster.
        cluster_number : int
            The ``Cluster`` identification number.  Distinguishes ``Cluster``
            from other ``Clusters`` in its ``ClusterGroup``

        """
        HyPath.__init__(self, clusterdata, pathdata, datetime,
                        clusterheader)
        # Initializes self.trajectories, self.trajcount, and self.directory
        HyGroup.__init__(self, trajectories)

        self.start_longitude = self.trajectories[0].loc[0, 'geometry'].x
        self.clusternumber = cluster_number

    def addgroups(self, other):
        """
        Prints notice before calling ``HyGroup.addgroups()``

        Parameters
        ----------
        other : ``TrajectoryGroup``
            Another ``TrajectoryGroup`` or ``Cluster``

        """

        print("Basic TrajectoryGroup created, cluster methods unavailable")

        # Initializes self.trajectories, self.trajcount, and self.directory
        new_tg = TrajectoryGroup(HyGroup.addgroups(self, other))

        return new_tg

    def calculate_vector(self):
        """
        Calculate mean bearing of ``Cluster`` path and bearings between
            timesteps

        """

        HyPath.calculate_vector(self)

    def calculate_distance(self):
        """
        Calculate the distance between timesteps fo the ``Cluster`` path and
            the cumulative distance at each time step

        """

        HyPath.calculate_distance(self)

    def pop(self, ind=-1, trajid=None):
        """
        Intercept ``HyGroup.pop()``.

        If a list of ``Trajectory`` instances is returned from
        ``HyGroup.pop()``, return a new ``TrajectoryGroup``,
        NOT a new ``Cluster``.

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


class ClusterGroup(object):
    """
    Contains all the ``Cluster``s produced in one ``HYSPLIT`` cluster analysis.

    """

    def __init__(self, clusters):
        """
        Initialize ``ClusterGroup`` object.

        Parameters
        ----------
        clusters : list of ``Cluster`` instances
            ``Cluster`` instances from the same HYSPLIT clustering run.
        """

        self.clusters = clusters

        self.clustercount = len(clusters)

        self.trajcount = sum([c.trajcount for c in self.clusters])

    def __getitem__(self, index):
        """
        Index or slice ``self.clusters`` to get a ``Cluster`` or
        ``ClusterGroup``, respectively.

        """

        newthing = self.clusters[index]

        try:
            newthing = ClusterGroup(newthing)
        except:
            pass

        return newthing
