from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as clr
import traj_accessory as ta
from trajgroup import TrajectoryGroup
import mapmaker as mm


class Cluster(TrajectoryGroup):
    """
    A special :subclass: of ``TrajectoryGroup`` for trajectories that have been
    clustered together using HYSPLIT's clustering process.

    Contains ``TrajectoryGroup`` attributes and functions, but also has
    ``Trajectory``-like attributes and functions associated with it, since
    a ``Cluster`` may be represented as a mean trajectory.
    """

    def __init__(self, traj_object_list, cluster_number, latitude, longitude):
        """
        Initialize ``Cluster`` object.

        Parameters
        ----------
        traj_object_list : list of ``Trajectory`` objects
            Trajectories that belong in the cluster.
        cluster_number : int
            The ``Cluster`` identification number.  Distinguishes ``Cluster``
            from other Clusters in its ``ClusterGroup``

        """
        # Initializes self.trajectories, self.trajcount, and self.directory
        TrajectoryGroup.__init__(self, traj_object_list)

        self.start_longitude = self.trajectories[0].longitude[0]
        self.clusternumber = cluster_number
        self.latitude = latitude
        self.longitude = longitude

        try:
            self.set_meanvar()
        except:
            print 'Unable to initialize mean variables'

    def __add__(self, other):
        """
        Prints notice before calling ``TrajectoryGroup.__add__()``

        Parameters
        ----------
        other : TrajectoryGroup

        """

        print "Basic TrajectoryGroup created, cluster methods unavailable"

        # Initializes self.trajectories, self.trajcount, and self.directory
        new_tg = TrajectoryGroup.__add__(self, other)

        return new_tg

    def set_vector(self):
        """
        Calculate mean bearing of ``Cluster`` path and bearings between
            timesteps

        """

        # Coordinates must be loaded first
        if not hasattr(self, 'latitude'):
            raise ValueError('self.latitude and self.longitude missing \n' +
                             'perform set_coordinates() first')

        self.meanvector, self.bearings = ta.tracemean_vector(self.latitude,
                                                             self.longitude)

    def set_distance(self):
        """
        Calculate the distance between timesteps fo the ``Cluster`` path and
            the cumulative distance at each time step

        """

        # Coordinates must be loaded first
        if not hasattr(self, 'latitude'):
            raise ValueError('self.latitude and self.longitude missing \n' +
                             'perform set_coordinates() first')
        self.distance = ta.distance_overearth(self.latitude, self.longitude)
        self.total_distance = ta.sum_distance(self.distance)

    def set_meanvar(self):
        """
        Calculate means and sums of several moisture-related variables.

        Means include the mean of every value along every trajectory in the
        Cluster and the mean of every trajectory's latest value.  Sums are
        the sum of every trajectory's latest value.

        Attributes set:
            mean_mf
            mean_w
            mean_q
            mean_raint0
            total_raint0
            mean_mft1
            total_mft1

        """

        mf_list = []
        q_list = []
        w_list = []
        raint0_list = []
        mft1_list = []

        for traj in self.trajectories:
            # Get trajectory attributes, if necessary setting them first
            if not hasattr(traj, 'specific_humidity'):
                traj.set_specifichumidity()

            if not hasattr(traj, 'mixing_ratio'):
                traj.set_mixingratio()

            if not hasattr(traj, 'moistureflux'):
                traj.calculate_moistureflux()

            mf_list.extend(traj.moistureflux[:-1].tolist())
            w_list.extend(traj.mixing_ratio.tolist())
            q_list.extend(traj.specific_humidity.tolist())
            raint0_list.append(traj.rainfall[0])
            mft1_list.append(traj.moistureflux[0])

        self.mean_mf = np.mean(mf_list)
        self.mean_w = np.mean(w_list)
        self.mean_q = np.mean(q_list)

        self.mean_raint0 = np.mean(raint0_list)
        self.total_raint0 = np.sum(raint0_list)
        self.mean_mft1 = np.mean(mft1_list)
        self.total_mft1 = np.sum(mft1_list)
        self.total_mf = np.sum(mf_list)

    def map_cluster_path(self, cavemap, color, lw=2, **kwargs):
        """
        Draw the map! DRAW IT NOW!

        Parameters
        ----------
        cavemap :  ``matplotlib`` ``Basemap`` instance
            Any basemap
        color : string or tuple of floats
            Any ``matplotlib``-accepted color
        lw : int
            Default 2.  The linewidth of the cluster path
        **kwargs
            Passed to ``mapmaker.traj_path()``.  Any ``Axes.plot()`` kwargs

        """

        cavemap = mm.traj_path(cavemap, self.longitude, self.latitude,
                               color, lw, **kwargs)


class ClusterGroup(object):
    """
    Class for processing and plotting Clusters
    """

    def __init__(self, cluster_object_list):
        """
        Initialize ClusterGroup object.

        Parameters
        ----------
        cluster_object_list : list of cluster objects
            Cluster objects from the same HYSPLIT clustering run.
        """

        self.totalclusters = len(cluster_object_list)
        self.clusters = cluster_object_list
        totaltraj = []
        for cluster in self.clusters:
            totaltraj.append(cluster.trajcount)
        self.totaltrajcount = sum(totaltraj)

    def map_clusters(self, cavemap, color_var='mean_mf', colors=None,
                     colormap=plt.cm.Blues, width_var='relative count',
                     lw=2, cnormalize=None, wnormalize=None, vmin=None,
                     vmax=None, **kwargs):
        """
        Plots mean cluster paths with color and width scaled to some variable.

        Creates a map where each cluster trajectory is plotted with a color
        determined by the mean variable value and a width by the trajectory
        count.  Color and width scaled according to user preferences.

        Parameters
        ----------
        cavemap : Basemap instance
            Initialize a basemap first using ``MapDesign.make_basemap()``
        color_var : string
            Default 'mean_mf'.  ['total_raint0'|'total_mft1'|'mean_q'|
            'mean_w'|'mean_mf'|'mean_mft1'|'mean_raint0'|'total_mf']
        color_min : int or float
            Default None.  The minimum value for color mapping.  If None,
            color_min will be the minimum value of the data.
        color_max : int or float
            Default None.  The maximum value for color mapping.  If None,
            color_max will be the maximum value of the data.
        color_rescale : string
            Default None.  ['sqrt'|ln'|'log']
            Determines how color data is transformed, if at all
        width_var : string
            Default 'count'.  ['absolute_count'|'relative_count'|
            'total_raint0'|'total_mft1'|'mean_q'|
            'mean_w'|'mean_mf'|'mean_mft1'|'mean_raint0'|'total_mf']
        width_rescale : string
            Default None.  [sqrt'|'ln'|'log']
            Determines how width data is transformed, if at all
        width_adjust : float
            Default 1.0.  Used to adjust linewidths to a reasonable value
            while maintaining relationships between linesizes.
        colormap : string
            Default 'blues'.  ['jet'|'blues'|'anomaly'|'heat'|'earth']
            Passed to a dictionary which retrieves the indicated colormap

        Returns
        -------
        colors : list of RGBA tuples
            The color of each cluster
        plot_order : list of ints
            The order of `colors` relative to ``self.clusters``

        """
        transform_dict = {'sqrt' : np.sqrt,
                          'log'  : np.log10,
                          'ln'   : np.ln}

        if width_var is None:
            wdata = [lw] * self.totalclusters
            plot_order = range(0, self.totalclusters)
        else:
            wdata = []
            for clus in self.clusters:
                if 'count' in width_var:
                    wdata.append(clus.trajcount)

                else:
                    wdata.append(getattr(clus, width_var))

            wdata = np.asarray(wdata, dtype=np.float64)

            if width_var is 'relative count':
                wdata = (wdata / self.totaltrajcount) * 100

            if wnormalize is not None:
                wdata = transform_dict[wnormalize](wdata)
            plot_order = np.argsort(wdata, kind='mergesort')

        if color_var is None:
            if colors is None:
                colors = mm.random_colors(self.totalclusters)
            else:
                try:
                    if len(colors) != self.totalclusters:
                        colors = [colors[0]] * self.totalclusters
                        print 'Number of colors and clusters do not match'
                except:
                    colors = [colors[0]] * self.totalclusters
        else:
            cdata = []
            for clus in self.clusteres:
                cdata.append(getattr(clus, color_var))

            cdata = np.asarray(cdata, dtype=np.float64)

            if cnormalize is not None:
                cdata = transform_dict[cnormalize](cdata)

            if vmin is not None:
                vmin = np.min(cdata)
            if vmax is not None:
                vmax = np.max(cdata)

            cnorm = clr.Normalize(vmin=vmin, vmax=vmax)
            scalarmap = cm.ScalarMappable(norm=cnorm, cmap=colormap)

            colors = []
            for c in cdata:
                color = scalarmap.to_rgba(c)
                colors.append(color)

        for i in plot_order:
            self.clusters[i].map_cluster_path(cavemap, colors[i], lw=wdata[i],
                                              **kwargs)

        return colors, plot_order

    def trajplot_clusterplot(self, trajmap, clusmap, colors=None, traj_lw=3,
                             cluster_lw='relative', **kwargs):

        """
        Plots trajectories on one map and cluster paths on another.

        Trajectories belonging to a cluster will be plotted as lines of the
        same color, and in the cluster plot the cluster mean path will also be
        the same color.  The linewidth of the cluster mean path will reflect
        the number of member trajectories unless specified otherwise.

        Parameters
        ----------
        trajmap : Basemap instance
            Initialize a basemap first using ``MapDesign.make_basemap()``
            Will have a plot of all trajectories
        clusmap : Basemap instance
            Initialize a basemap first using ``MapDesign.make_basemap()``
            Will have a plot of all clusters
        colors : list of strings or RGB tuples
            Default None.  Length of list should match self.totalclusters.
            If None, a list of random colors will be generated.
        cluster_lw : int or string
            Default 'relative'.  [int|'relative'|'absolute'].
            The linewidth of the cluster mean path is set with a constant
            value, or is set with a value that depends on the relative or
            absolute number of member trajectories in the cluster
        traj_lw : int
            Default 3.  The linewdith of each trajectory.

        Returns
        -------
        colors : list of RGBA tuples
            The color used for each cluster.

        """

        # Set colors if none provided
        if colors is None:
            mm.random_colors(self.totalclusters)

        # Set cluster mean path width by member trajectory count if not given
        clus_lw = []

        if cluster_lw is 'relative':
            for cluster in self.clusters:
                clus_lw.append((cluster.trajcount /
                                float(self.totaltrajcount)) * 100.0)
        elif cluster_lw is 'absolute':
            for cluster in self.clusters:
                clus_lw.append(cluster.trajcount)
        else:
            clus_lw = [cluster_lw]
            clus_lw = cluster_lw * self.totalclusters

        # Plot
        for cluster, color, lw in zip(self.clusters, colors, cluster_lw):
            cluster.map_cluster_path(clusmap, color, lw=lw, **kwargs)

            for traj in cluster.trajectories:
                traj.map_traj_path(trajmap, color=color, lw=traj_lw, **kwargs)

        return colors
