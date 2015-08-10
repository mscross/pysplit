from __future__ import division, print_function
from trajgroup import TrajectoryGroup
from hypath import HyPath
from hygroup import HyGroup


class Cluster(HyPath, HyGroup):
    """
    A special :subclass: of ``HyGroup`` for trajectories that have been
    clustered together using HYSPLIT's clustering process.

    Contains ``HyGroup`` attributes and functions, but also has
    ``Trajectory``-like attributes and functions associated with it, since
    a ``Cluster`` may be represented as a mean trajectory.
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

        self.start_longitude = self.trajectories[0].longitude[0]
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

        HyPath.set_vector(self)

    def set_distance(self):
        """
        Calculate the distance between timesteps fo the ``Cluster`` path and
            the cumulative distance at each time step

        """

        HyPath.set_distance(self)


class ClusterGroup(object):
    """
    Class for processing and plotting member ``Cluster`` instances
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

        self.trajcount = sum([len(c.trajectories) for c in self.clusters])

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
