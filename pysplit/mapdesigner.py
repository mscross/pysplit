from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from .maplabeller import map_labeller, labelfile_reader


class MapDesign(object):
    """
    Class for holding map design elements

    """

    def __init__(self, mapcorners, standard_pm, projection='cyl',
                 mapcolor='light', maplabels=None, area_threshold=10000,
                 resolution='c', zborder=14, lat_labels=['left'],
                 lon_labels=['top'], latlon_labelspacing=(10, 20),
                 latlon_fs=20, latlon_spacing=(10, 20), drawstates=False,
                 drawoutlines=True, draw_latlons=True):
        """
        Initialize ``MapDesign`` instance.

        Parameters
        ----------
        mapcorners : list of floats
            Used to construct the map view for 'conic' and 'cyl' projections.
            Lower left longitude, latitude; upper right longitude, latitude.
        standard_pm : list of floats
            For cylindrical and conic projections, the list creates standard
                parallels and meridians
                (``lon_0``, ``lat_0``, ``lat_1``, ``lat_2``).
            For orthographic projection, ``lon_0`` and ``lat_0`` only
                are required.  Sets the view from above the earth.
            For polar projections, ``lon_0`` indicates the longitude that
                will be oriented N-S. ``lat_0`` is replaced by
                the ``boundinglat``, the lowest latitude that should appear
                on the map.  ``lat_1`` and ``lat_2`` not required
        projection : string
            Indicates map projection.  Default 'cyl'.
                'cyl' : Equidistant cylindrical
                'cea' : Equal Area cylindrical
                'lcc' : Lambert Conformal Conic
                'aea' : Albers Equal Area Conic
                'ortho' : Orthographic (globe)
                'npstere' : North polar steroegraphic (conformal)
                'spstere' : South polar steroegraphic (conformal)
                'nplaea' : North polar azimuthal (equal area)
                'splaea' : South polar azimuthal (equal area)
        mapcolor : string
            Default 'light'. The map grayscheme.
            ['light'|'medium'|'dark'|None]
            Not available for 'ortho' projections
        maplabels : tuple of strings
            Default ``None``.
            (Label group, label file full/relative path, optional: zorder).
            Label group is a list of any or all of:
            ['sea', 'city', 'country','ocean','place']
            Label zorder defaults to 15.
        area_threshold : int
            Default 10000.  The minimum surface area a feature must have to
            be drawn on the map.
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
        latlon_labelspacing : int or float
            Default (10, 20).  Degrees between (latitude, longitude) labels
        latlon_fs : int or float
            Default 20.  Font size of latitude, longitude labels.
        latlonspacing : int or float
            Default (10, 20).  Degrees between plotted lines of latitude.
        drawstates : Boolean
            Default False.  Draw state outlines on ``Basemap``.
        drawoutlines : Boolean
            Default True.  Draw country and coastal outlines on ``Basemap``.
        draw_latlons : Boolean
            Default True.  Draw and label lines of latitude and longitude.

        """

        # Initialize
        self.mapcorners = mapcorners
        self.standard_pm = standard_pm

        self.mapcolor = mapcolor
        self.coloropts = {'light': {'water': None,
                                    'land': '0.95'},
                          'medium': {'water': '0.625',
                                     'land': '0.775'},
                          'dark': {'water': '0.3',
                                   'land': '0.75'}}

        self.area_threshold = area_threshold
        self.resolution = resolution
        self.zborder = zborder

        # Initialize projection
        self._set_projection(projection)

        self._set_latlonlabels(lat_labels, lon_labels)

        self.latlon_fs = latlon_fs

        self.latspacing = latlon_spacing[0]
        self.lonspacing = latlon_spacing[1]
        self.latstep = latlon_labelspacing[0]
        self.lonstep = latlon_labelspacing[1]

        self.drawstates = drawstates
        self.drawoutlines = drawoutlines
        self.draw_latlons = draw_latlons

        # Try to set label attributes
        if maplabels is not None:

            self.labels, self.labelstyle = labelfile_reader(maplabels[1])
            self.labelgroup = maplabels[0]

            # Label zorder optional, default 15 if none given.
            try:
                self.label_zorder = maplabels[2]
            except:
                self.label_zorder = 15
        else:
            self.labels = None

    def _set_latlonlabels(self, lat_labels, lon_labels):
        """
        Edit latitude, longitude labelling preferences.

        Parameters
        ----------
        lat_labels : list of strings
            The sides of the map with longitude labels.
        lon_labels : list of strings
            The sides of the map with latitude labels.

        """

        meridian_labels = [0, 0, 0, 0]
        parallel_labels = [0, 0, 0, 0]

        ind_dict = {'left': 0,
                    'right': 1,
                    'top': 2,
                    'bottom': 3}

        for la in lat_labels:
            parallel_labels[ind_dict[la]] = 1
        self.parallel_labels = parallel_labels

        for lo in lon_labels:
            meridian_labels[ind_dict[lo]] = 1
        self.meridian_labels = meridian_labels

    def _set_projection(self, projection):
        """
        Set the projection.  Defaults to 'cyl'.

        Parameters
        ----------
        projection : string
        Indicates which projection to use.  Default 'cyl'
            'cyl' : Equidistant cylindrical
            'cea' : Equal Area cylindrical
            'lcc' : Lambert Conformal Conic
            'aea' : Albers Equal Area Conic
            'ortho' : Orthographic (globe)
            'npstere' : North polar steroegraphic (conformal)
            'spstere' : South polar steroegraphic (conformal)
            'nplaea' : North polar azimuthal (equal area)
            'splaea' : South polar azimuthal (equal area)

        """

        available_proj = {'cyl': 'Equidistant cylindrical',
                          'cea': 'Equal Area cylindrical',
                          'lcc': 'Lambert Conformal Conic',
                          'aea': 'Albers Equal Area',
                          'ortho': 'Orthographic',
                          'npstere': 'North Polar Stereographic',
                          'spstere': 'South polar steroegraphic',
                          'nplaea': 'North polar azimuthal',
                          'splaea': 'South polar azimuthal'}

        if projection in available_proj:
            self.projection = projection
        else:
            self.projection = 'cyl'
            print('Projection not recognized, defaulting to `cyl`.')

    def make_basemap(self, ax=None, figsize=(10, 10)):
        """
        Takes the MapDesign attributes plus a figure size and creates a map
        on which data can be plotted.

        Parameters
        ----------
        ax : axes instance
            Default None, figure and axis will be created.  Otherwise,
            basemap will be created on given axis.
        figsize : tuple of ints
            Default (10, 10). The size of the figure in inches.  Only
            used if ``ax`` is ``None``.

        Returns
        -------
        basemap : ``Basemap`` instance
            A map ready for data plotting.  Can access axis and figure
            via ``basemap.ax`` and ``basemap.ax.get_figure()``, respectively.

        """

        if ax is None:
            # Create figure instance
            fig, ax = plt.subplots(1, 1, figsize=figsize)

        # Labels are left, right, top, bottom
        meridian_labels = self.meridian_labels
        parallel_labels = self.parallel_labels

        if self.projection is 'lcc':
            # Lambert conformal conic
            basemap = Basemap(llcrnrlon=self.mapcorners[0],
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

        elif self.projection is 'aea':
            # Albers equal area conic
            basemap = Basemap(llcrnrlon=self.mapcorners[0],
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

        elif self.projection is 'cea':
            # equal area cylindrical
            basemap = Basemap(llcrnrlon=self.mapcorners[0],
                              llcrnrlat=self.mapcorners[1],
                              urcrnrlon=self.mapcorners[2],
                              urcrnrlat=self.mapcorners[3],
                              projection='cea',
                              area_thresh=self.area_threshold,
                              resolution=self.resolution,
                              ax=ax)

        elif self.projection is 'ortho':
            # the earth
            basemap = Basemap(projection='ortho',
                              lon_0=self.standard_pm[0],
                              lat_0=self.standard_pm[1],
                              area_thresh=self.area_threshold,
                              resolution=self.resolution,
                              ax=ax)

            meridian_labels = [0, 0, 0, 0]
            parallel_labels = [0, 0, 0, 0]

        elif self.projection[1:] == 'plaea' or self.projection[1:] == 'pstere':
            # Polar steroegraphic (conformal) or polar azimuthal(equal-area)
            basemap = Basemap(projection=self.projection,
                              boundinglat=self.standard_pm[1],
                              lon_0=self.standard_pm[0],
                              area_thresh=self.area_threshold,
                              resolution=self.resolution,
                              ax=ax)

        else:
            # Projection is 'cyl', Basemap's default equidist. cyl projection
            basemap = Basemap(llcrnrlon=self.mapcorners[0],
                              llcrnrlat=self.mapcorners[1],
                              urcrnrlon=self.mapcorners[2],
                              urcrnrlat=self.mapcorners[3],
                              area_thresh=self.area_threshold,
                              resolution=self.resolution,
                              ax=ax)

        # Draw labels
        if self.labels is not None:
            map_labeller(basemap, self.labelgroup,
                         self.labels, self.labelstyle,
                         self.label_zorder)

        # Draw countries, states, coastlines, and map boundary
        if self.drawstates:
            basemap.drawstates(zorder=14)
        if self.drawoutlines:
            basemap.drawcountries(zorder=self.zborder)
            basemap.drawcoastlines(zorder=self.zborder)

        if self.draw_latlons:
            lat = -90
            if self.latstep == 20:
                lat = -80

            # Draw and label lines of longitude
            if self.lonstep > self.lonspacing:

                basemap.drawmeridians(np.arange(-180, 180, self.lonspacing),
                                      labels=[0, 0, 0, 0], zorder=11)
                basemap.drawmeridians(np.arange(-180, 180, self.lonstep),
                                      labels=meridian_labels, zorder=11,
                                      fontsize=self.latlon_fs)
            else:
                basemap.drawmeridians(np.arange(-180, 180, self.lonstep),
                                      labels=meridian_labels, zorder=11,
                                      fontsize=self.latlon_fs)

            # Draw and label lines of latitude
            if self.latstep > self.latspacing:

                basemap.drawparallels(np.arange(-90, 90, self.latspacing),
                                      labels=[0, 0, 0, 0], zorder=11)
                basemap.drawparallels(np.arange(lat, 90, self.latstep),
                                      labels=parallel_labels, zorder=11,
                                      fontsize=self.latlon_fs)
            else:
                basemap.drawparallels(np.arange(lat, 90, self.latstep),
                                      labels=parallel_labels, zorder=11,
                                      fontsize=self.latlon_fs)

        # Map color defaults to white for ortho projection
        if self.mapcolor is not None:
            if self.projection is not 'ortho':
                colors = self.coloropts[self.mapcolor]

                basemap.drawmapboundary(zorder=16, fill_color=colors['water'])
                basemap.fillcontinents(color=colors['land'], zorder=12,
                                       alpha=0.85, lake_color=colors['water'])
            else:
                basemap.fillcontinents(color='0.99', zorder=12, alpha=0.85)

        return basemap
