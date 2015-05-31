from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import maplabeller as ml


class MapDesign(object):
    """
    Class for holding map design elements

    """

    def __init__(self, mapcorners, standard_pm, projection='cyl',
                 mapcolor='light', shapefiles=[], maplabels=None,
                 area_threshold=10000, resolution='c', zborder=14,
                 lat_labels=['left'], lon_labels=['top'], lat_labelspacing=10,
                 lon_labelspacing=20, latlon_fs=20, latspacing=10,
                 lonspacing=20):
        """
        Initialize ``MapDesign`` instance.

        Parameters
        ----------
        mapcorners : list of floats
            Used to construct the map view for 'conic' and 'cyl' projections.
            Lower left longitude, latitude; upper right longitude, latitude.
        standard_pm : list of floats
            For `cyl` and `conic` projections, the list creates standard
                parallels and meridians
                (``lon_0``, ``lat_0``, ``lat_1``, ``lat_2``).
            For 'orthographic' projection, ``lon_0`` and ``lat_0`` only
                are required.  Sets the view from above the earth.
            For polar projections, ``lon_0`` indicates the longitude that
                will be oriented N-S. ``lat_0`` is replaced by
                the ``boundinglat``, the lowest latitude that should appear
                on the map.  ``lat_1`` and ``lat_2`` not required
        projection : string
            Indicates which projection to use.  Default 'cyl'
                'cyl' : Equidistant cylindrical
                'eacyl' : Equal Area cylindrical
                'cconic' : Lambert Conformal Conic
                'eaconic' : Albers Equal Area Conic
                'ortho' : Orthographic (globe)
                'npstere' : North polar steroegraphic (conformal)
                'spstere' : South polar steroegraphic (conformal)
                'nplaea' : North polar azimuthal (equal area)
                'splaea' : South polar azimuthal (equal area)
        mapcolor : string
            Default 'light'. The map grayscheme. ['light'|'medium'|'dark']
            Not available for 'ortho' projections
        shapefiles : list of tuples of strings
            Default is [].  (File, color, linewidth)
            r'C:\programming\shapefiles\New_Shapefile'
        maplabels : tuple of strings
            Default ``None``.
            (Label group, label file full/relative path, optional: zorder).
            Determines what label groups are applied, if any:
            Label group choices:  ['all_1'|'all_2'|important'|
                'justcave'|'cave']
            If path to label file does not exist, the user will be presented
            with several options, including one to make and edit a label file.
            Label zorder defaults to 15.
        area_threshold : int
            Default 10000.  The minimum surface area a feature must have to
            be plotted on the map.
        resolution : char
            Default 'c'.  ['c'|'l'|'i'|'h'|'f'].
            Crude, low, intermediate, high, full. The relative resolution of
            map boundaries.  Drops off by about 80 percent between datasets.
        zborder : int
            Default 14. The zorder of country and coastal outlines
        lat_labels : list of strings
            Default ['right'].
            The sides of the map that should have the latitudes labelled.
        lon_labels : list of strings
            Default ['top'].
            The sides of the map that should have the longitudes labelled.
        lat_labelspacing : int or float
            Default 10.  Degrees between latitude labels
        lon_labelspacing : int or float
            Default 20.  Degrees between longitude labels.
        latlon_fs : int or float
            Default 20.  Font size of latitude, longitude labels.
        latspacing : int or float
            Default 10.  Degrees between plotted lines of latitude.
        lonspacing : int or float
            Default 20.  Degrees between plotted lines of longitude.

        """

        # Initialize
        self.mapcorners = mapcorners
        self.standard_pm = standard_pm

        self.mapcolor = mapcolor
        self.coloropts = {'light' : {'water' : None,
                                     'land' : '0.95'},
                          'medium' : {'water' : '0.625',
                                      'land' : '0.775'},
                          'dark' : {'water' : '0.3',
                                    'land' : '0.75'}}

        self.shapefiles = shapefiles
        self.area_threshold = area_threshold
        self.resolution = resolution
        self.zborder = zborder

        # Initialize projection
        self.edit_projection(projection)

        self.edit_latlonlabels(lat_labels=lat_labels, lon_labels=lon_labels,
                               lat_labelspacing=lat_labelspacing,
                               lon_labelspacing=lon_labelspacing,
                               latlon_fs=latlon_fs)

        self.latspacing = latspacing
        self.lonspacing = lonspacing

        # Try to set label attributes
        if maplabels is not None:
            self.labels, self.labelstyle = ml.labelfile_reader(maplabels[1])
            self.labelgroup = maplabels[0]

            # Label zorder optional, default 15 if none given.
            try:
                self.label_zorder = maplabels[2]
            except:
                self.label_zorder = 15
        else:
            self.clear_labels()

    def view_prefs(self):
        """
        Create a table of the current map design elements (i.e. attributes)
        and their values.

        """

        pref_list = self.__dict__.keys()

        for pref, num in zip(pref_list, range(1, len(pref_list) + 1)):
            print '\t', num, '. ', pref, ' : ', getattr(self, pref)

        print '\n'

    def edit_latlonlabels(self, lat_labels=None, lon_labels=None,
                          lat_labelspacing=None, lon_labelspacing=None,
                          latlon_fs=None):
        """
        Edit latitude, longitude labelling preferences.

        Parameters
        ----------
        lat_labels : list of strings
            The sides of the map that should have the latitudes labelled.
        lon_labels : list of strings
            The sides of the map that should have the longitudes labelled.
        lat_labelspacing : int or float
            Degrees between latitude labels
        lon_labelspacing : int or float
            Degrees between longitude labels.
        latlon_fs : int or float
            Font size of latitude, longitude labels.

        """

        meridian_labels = [0, 0, 0, 0]
        parallel_labels = [0, 0, 0, 0]

        ind_dict = {'left' : 0,
                    'right' : 1,
                    'top' : 2,
                    'bottom' : 3}

        if lat_labels is not None:
            if 'none' not in lat_labels:
                for la in lat_labels:
                    parallel_labels[ind_dict[la]] = 1
            self.parallel_labels = parallel_labels

        if lon_labels is not None:
            if 'none' not in lon_labels:
                for lo in lon_labels:
                    meridian_labels[ind_dict[lo]] = 1
            self.meridian_labels = meridian_labels

        if lat_labelspacing is not None:
            self.latstep = lat_labelspacing

        if lon_labelspacing is not None:
            self.lonstep = lon_labelspacing

        if latlon_fs is not None:
            self.latlon_fs = latlon_fs

    def edit_latlonspacing(self, lonspacing=None, latspacing=None):
        """
        Change the spacing between the lines of latitude, longitude drawn on
        the map.

        Parameters
        ----------
        lonspacing : int or float
            Default 20.  Degrees between plotted lines of longitude.
        latspacing : int or float
            Default 10.  Degrees between plotted lines of latitude.

        """

        if lonspacing is not None:
            self.lonspacing = lonspacing

        if latspacing is not None:
            self.latspacing = latspacing

    def edit_resolution(self, resolution=None, area_threshold=None):
        """
        Adjust the map resolution and area threshold for plotting.
        Attributes will only be adjusted if not ``None``.

        Parameters
        ----------
        resolution : string
            Default 'c'.  ['c'|'l'|'i'|'h'|'f'].
            Crude, low, intermediate, high, full. The relative resolution of
            map boundaries.  Drops off by about 80 percent between datasets.
        area_threshold : int or string
            Default ``None``.  ['auto'|1|10|100|1000|10000]
            The minimum surface area a feature must have to
            be plotted on the map.  If 'auto', an area threshold
            will be selected based on the current map resolution.

        """

        if resolution is not None:
            self.resolution = resolution

        if area_threshold is not None:
            if area_threshold == 'auto':
                self.area_threshold = None
            else:
                self.area_threshold = area_threshold

    def edit_mapcorners(self, mapcorners):
        """
        Update the map extent.

        Parameters
        ----------
        mapcorners : list of floats
            Used to construct the map view for 'conic' and 'cyl' projections.
            Lower left longitude, latitude; upper right longitude, latitude.

        """

        self.mapcorners = mapcorners

    def edit_standard_pm(self, standard_pm):
        """
        Update the standard parallels and meridians.

        Parameters
        ----------
        standard_pm : list of floats
            For `cyl` and `conic` projections, the list creates standard
                parallels and meridians
                (``lon_0``, ``lat_0``, ``lat_1``, ``lat_2``).
            For 'orthographic' projection, ``lon_0`` and ``lat_0`` only
                are required.  Sets the view from above the earth.
            For polar projections, ``lon_0`` indicates the longitude that
                will be oriented N-S. ``lat_0`` is replaced by
                the ``boundinglat``, the lowest latitude that should appear
                on the map.  ``lat_1`` and ``lat_2`` not required

        """

        self.standard_pm = standard_pm

    def set_mapcolor(self, mapcolor):
        """
        Change the map grayscheme.

        Parameters
        ----------
        mapcolor : string
            Grayscheme option.  ['light'|'medium'|'dark']

        """

        options = ['light', 'medium', 'dark']

        if mapcolor not in options:
            print 'Unrecognized mapcolor ', mapcolor
            self.mapcolor = 'light'
        else:
            self.mapcolor = mapcolor

    def edit_shapefiles(self, shp=None, delete_shp=None, overwrite_shps=False):
        """
        Edit the list of shapefiles.

        Parameters
        ----------
        shp : tuple of list of tuples of strings
            Default ``None``.  Tuples  are (File, color, linewidth)
            r'C:\programming\shapefiles\New_Shapefile'.
            Used to extend or overwrite current list of shapefiles.
        delete_shp : int or list of ints
            Default ``None``.  The indices of shapefiles to delete.
        overwrite_shps : Boolean
            Default ``False``.  Indicates whether to extend the current
            shapefiles attributes with new shapes or to overwrite.

        """

        # Remove shapefile(s) by index
        if delete_shp is not None:
            try:
                for ds in delete_shp:
                    del self.shapefiles[ds]
            except:
                del self.shapefiles[delete_shp]

        # Make list of shapes the new shapefile list or add it to  current list
        if shp is not None:
            if type(shp) is not list:
                shp = [shp]
            if overwrite_shps:
                self.shapefiles = shp
            else:
                self.shapefiles.extend(shp)
        else:
            if overwrite_shps:
                self.shapefiles = []

    def edit_projection(self, projection):
        """
        Update the projection.  Defaults to 'cyl'.

        Parameters
        ----------
        projection : string
        Indicates which projection to use.  Default 'cyl'
            'cyl' : Equidistant cylindrical
            'eacyl' : Equal Area cylindrical
            'cconic' : Lambert Conformal Conic
            'eaconic' : Albers Equal Area Conic
            'ortho' : Orthographic (globe)
            'npstere' : North polar steroegraphic (conformal)
            'spstere' : South polar steroegraphic (conformal)
            'nplaea' : North polar azimuthal (equal area)
            'splaea' : South polar azimuthal (equal area)

        """

        available_proj = {'cyl' : 'Equidistant cylindrical',
                          'eacyl' : 'Equal Area cylindrical',
                          'cconic' : 'Lambert Conformal Conic',
                          'eaconic' : 'Albers Equal Area',
                          'ortho' : 'Orthographic',
                          'npstere' : 'North Polar Stereographic',
                          'spstere' : 'South polar steroegraphic',
                          'nplaea' : 'North polar azimuthal',
                          'splaea' : 'South polar azimuthal'}

        if projection in available_proj:
            self.projection = projection
        else:
            self.projection = 'cyl'
            # Error message
            print 'Projection not recognized.'
            print 'Choose a projection from the left column:'
            for proj in available_proj:
                print '\t', proj, '\t', available_proj[proj]

    def edit_labels(self, labelpath=None, labelgroup=None, label_zorder=None):
        """
        Change which labels are applied to the map.  Attributes will change
        only if not None.

        Parameters
        ----------
        labelpath : string
            Default ``None``.  Full or relative path to labelfile location
        labelgroup : string
            Default ``None``.  ['all_1'|'all_2'|important'|'justcave'|'cave']
            The set of labels to apply.  'all_1' has caves instead of cities,
            'all_2' has cities instead of caves.

        """

        if labelpath is not None:
            self.labels, self.labelstyle = ml.labelfile_reader(labelpath)

        if labelgroup is not None:
            self.labelgroup = labelgroup

        if label_zorder is not None:
            self.label_zorder = label_zorder

    def clear_labels(self):
        """
        Reset label path and group information to None, label zorder to 15.

        """

        self.labels = None
        self.labelstyle = None
        self.labelgroup = None
        self.label_zorder = 15

    def edit_border_zorder(self, zborder):
        """
        Adjust the zorder of the country and coastal outlines.

        Parameters
        ----------
        zborder : int
            The zorder of country and coastal outlines.

        """

        self.zborder = zborder

    def make_basemap(self, ax=None, figsize=(15, 15)):
        """
        Takes the MapDesign attributes plus a figure size and creates a map
        on which data can be plotted.

        Parameters
        ----------
        ax : axes instance
            Default None, figure and axis will be created.  Otherwise,
            basemap will be created on given axis.
        figsize : tuple of ints
            Default (15, 15). The size of the figure in inches.  Only
            used if ax is ``None``.

        Returns
        -------
        cavemap : ``Basemap`` instance
            A map ready for data plotting.  Can access axis and figure
            via ``cavemap.ax`` and ``cavemap.ax.get_figure()``, respectively.

        """

        if ax is None:
            # Create figure instance
            fig, ax = plt.subplots(1, 1, figsize=figsize)

        # Labels are left, right, top, bottom
        meridian_labels = self.meridian_labels
        parallel_labels = self.parallel_labels

        if self.projection is 'cconic':
            # Lambert conformal conic
            cavemap = Basemap(llcrnrlon=self.mapcorners[0],
                              llcrnrlat=self.mapcorners[1],
                              urcrnrlon=self.mapcorners[2],
                              urcrnrlat=self.mapcorners[3],
                              projection='lcc',
                              lat_1=self.standard_pm[2],
                              lat_2=self.standard_pm[3],
                              lon_0=self.standard_pm[0],
                              area_thresh=self.area_threshold,
                              resolution=self.resolution,
                              ax=ax)

        elif self.projection is 'eaconic':
            # Albers equal area conic
            cavemap = Basemap(llcrnrlon=self.mapcorners[0],
                              llcrnrlat=self.mapcorners[1],
                              urcrnrlon=self.mapcorners[2],
                              urcrnrlat=self.mapcorners[3],
                              projection='aea',
                              lat_1=self.standard_pm[2],
                              lat_2=self.standard_pm[3],
                              lon_0=self.standard_pm[0],
                              lat_0=self.standard_pm[1],
                              area_thresh=self.area_threshold,
                              resolution=self.resolution,
                              ax=ax)

        elif self.projection is 'eacyl':
            # equal area cylindrical
            cavemap = Basemap(llcrnrlon=self.mapcorners[0],
                              llcrnrlat=self.mapcorners[1],
                              urcrnrlon=self.mapcorners[2],
                              urcrnrlat=self.mapcorners[3],
                              projection='cea',
                              area_thresh=self.area_threshold,
                              resolution=self.resolution,
                              ax=ax)

        elif self.projection is 'ortho':
            # the earth
            cavemap = Basemap(projection='ortho',
                              lon_0=self.standard_pm[0],
                              lat_0=self.standard_pm[1],
                              area_thresh=self.area_threshold,
                              resolution=self.resolution,
                              ax=ax)

            meridian_labels = [0, 0, 0, 0]
            parallel_labels = [0, 0, 0, 0]

        elif self.projection[1:] == 'plaea' or self.projection[1:] == 'pstere':
            # Polar steroegraphic (conformal) or polar azimuthal(equal-area)
            cavemap = Basemap(projection=self.projection,
                              boundinglat=self.standard_pm[1],
                              lon_0=self.standard_pm[0],
                              area_thresh=self.area_threshold,
                              resolution=self.resolution,
                              ax=ax)

        else:
            # Projection is 'cyl', Basemap's default equidist. cyl projection
            cavemap = Basemap(llcrnrlon=self.mapcorners[0],
                              llcrnrlat=self.mapcorners[1],
                              urcrnrlon=self.mapcorners[2],
                              urcrnrlat=self.mapcorners[3],
                              area_thresh=self.area_threshold,
                              resolution=self.resolution,
                              ax=ax)

        # Draw labels
        if self.labels is not None and self.labelgroup is not None:
                cavemap = ml.map_labeller(cavemap, ax, self.labelgroup,
                                          self.labels, self.label_zorder,
                                          self.labelstyle)

        if len(self.shapefiles) > 0:
            # plot shapefiles
            for shp in self.shapefiles:
                cavemap.readshapefile(shp[0], 'hella_shapes', zorder=15,
                                      color=shp[1], linewidth=shp[2])

        # Draw countries, states, coastlines, and map boundary
        cavemap.drawcountries(zorder=self.zborder)
        cavemap.drawstates(zorder=14)
        cavemap.drawcoastlines(zorder=self.zborder, linewidth=1.5)

        if self.latstep == 20:
            lat = -80
        else:
            lat = -90

        # Draw and label lines of longitude
        if self.lonstep > self.lonspacing:

            cavemap.drawmeridians(np.arange(-180, 180, self.lonspacing),
                                  labels=[0, 0, 0, 0], zorder=11)
            cavemap.drawmeridians(np.arange(-180, 180, self.lonstep),
                                  labels=meridian_labels, zorder=11,
                                  fontsize=self.latlon_fs)
        else:
            cavemap.drawmeridians(np.arange(-180, 180, self.lonstep),
                                  labels=meridian_labels, zorder=11,
                                  fontsize=self.latlon_fs)

        # Draw and label lines of latitude
        if self.latstep > self.latspacing:

            cavemap.drawparallels(np.arange(-90, 90, self.latspacing),
                                  labels=[0, 0, 0, 0], zorder=11)
            cavemap.drawparallels(np.arange(lat, 90, self.latstep),
                                  labels=parallel_labels, zorder=11,
                                  fontsize=self.latlon_fs)
        else:
            cavemap.drawparallels(np.arange(lat, 90, self.latstep),
                                  labels=parallel_labels, zorder=11,
                                  fontsize=self.latlon_fs)

        # Map color defaults to white for ortho projection
        if self.projection != 'ortho':
            colors = self.coloropts[self.mapcolor]

            cavemap.drawmapboundary(zorder=16, fill_color=colors['water'])
            cavemap.fillcontinents(color=colors['land'], zorder=12, alpha=0.85,
                                   lake_color=colors['water'])
        else:
            cavemap.fillcontinents(color='0.99', zorder=12, alpha=0.85)

        return cavemap
