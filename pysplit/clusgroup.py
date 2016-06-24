from __future__ import division, print_function
from trajgroup import TrajectoryGroup
from hypath import HyPath
from hygroup import HyGroup


def print_clusterprocedure():
    """Print clustering guide."""
    print("""
          In ``PySPLIT``
          1. Create ``TrajectoryGroup`` with desired set of trajectories
          2. ``TrajectoryGroup.make_infile()``

          In ``HYSPLIT``
          3. Trajectory --> Special Runs --> Clustering --> Standard
          4. Adjust clustering parameters, endpoints (trajectory) folder, and
             working folder (where output will be stored, where INFILE lives)
          5. ``Run cluster analysis``
          6. Determine and set appropriate number of clusters
          7. Assign trajectories to clusters (``Run``)
          8. ``Display Means``, ``Display Clusters``, ``Quit``

          In ``PySPLIT``
          9. ``spawn_clusters()``""")


class Cluster(HyPath, HyGroup):
    """
    A special :subclass: of both ``HyGroup`` and ``HyPath``.

    Clusters contain both trajectories and mean path information.  The mean
    path and the trajectory composition is determined by ``HySPLIT``.

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

        HyGroup.__init__(self, trajectories)

        self.start_longitude = self.trajectories[0].data.loc[0, 'geometry'].x
        self.clusternumber = cluster_number

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
        Add a ``HyGroup`` to this ``Cluster`` instance.

        Parameters
        ----------
        other : ``HyGroup``
            Another ``TrajectoryGroup`` or ``Cluster``.  May or may not
            contain some of the same ``Trajectory`` instances.

        Returns
        -------
        A new ``TrajectoryGroup`` containing the union of the sets
        of ``Trajectory`` instances.

        """
        return TrajectoryGroup(HyGroup.__add__(self, other))

    def __sub__(self, other):
        """
        Subtract a ``HyGroup`` from this ``Cluster`` instance.

        Parameters
        ----------
        other : ``HyGroup``
            Another ``Cluster`` or ``TrajectoryGroup``

        Returns
        -------
        A new ``TrajectoryGroup`` containing the set difference betwee
        the sets of ``Trajectory`` instances.

        """
        return TrajectoryGroup(HyGroup.__sub__(self, other))

    def calculate_vector(self):
        """
        Calculate bearings.

        Find mean bearing of ``Cluster`` path and bearings between
            timesteps.

        """
        HyPath.calculate_vector(self)

    def calculate_distance(self):
        """
        Calculate distances.

        Find the distance between timesteps of the ``Cluster`` path and
            the cumulative distance at each time step.

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
    Group of ``Cluster`` instances.

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
        Get ``Cluster`` or ``ClusterGroup``.

        Index or slice ``self.clusters`` to get a ``Cluster`` or
        ``ClusterGroup``, respectively.

        """
        newthing = self.clusters[index]

        try:
            newthing = ClusterGroup(newthing)
        except:
            pass

        return newthing
