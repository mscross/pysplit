import numpy as np
import math
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
from mpl_toolkits.basemap import Basemap
import maplabeller as ml

"""
To Do:
    -draw parallels and meridians twice
    -latspacing
    -lonspacing
    -if 20 is chosen as latlabelspacing, edit arange start -90 to -80
    -create medium background

"""


class MapDesign(object):
    """
    Class for holding map design elements

    """

    def __init__(self, mapcorners, standard_pm, projection='cyl',
                 mapcolor='light', shapefiles=[], maplabels=None,
                 area_threshold=10000, resolution='c', zborder=14,
                 lat_labels=['right'], lon_labels=['top'], lat_labelspacing=10,
                 lon_labelspacing=20, latlon_fs=20):
        """
        Initialize MapDesign object.

        Parameters
        ----------
        mapcorners : list of floats
            Used to construct the map view for conic and cyl projections.
            Lower left longitude, latitude; upper right longitude, latitude.
        standard_pm : list of floats
            For cyl and conic projections, the list creates standard parallels
                and meridians (lon_0, lat_0, lat_1, lat_2).
            For orthographic projection, lon_0 and lat_0 only are required.
                Sets the view from above the earth.
            For polar projections, lon_0 indicates the longitude that will be
                oriented N-S. lat_0 is replaced by the boundinglat,
                the lowest latitude that should appear on the map.
                lat_1 and lat_2 not required

        Keyword Arguments
        -----------------
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
            Default 'light' The map grayscheme. ['light'|'dark']
            'Dark' is not available for 'ortho' projections
        shapefiles : list of tuples of strings
            Default is [].  (File, color, linewidth)
            r'C:\programming\shapefiles\New_Shapefile'
        maplabels : tuple of strings
            Default None.
            (Label group, label file full/relative path, optional: zorder).
            Determines what label groups are applied, if any:
            Label group choices:  ['all_1'|'all_2'|important'|'justcave'|'cave']
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

        """

        # Initialize
        self.mapcorners = mapcorners
        self.standard_pm = standard_pm
        self.mapcolor = mapcolor
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

        # Try to set label attributes
        if maplabels is not None:
            self.labels = ml.labelfile_reader(maplabels[1])
            self.labelgroup = maplabels[0]

            # Label zorder optional, default 15 if none given.
            try:
                self.label_zorder = maplabels[2]
            except:
                self.label_zorder = 15
        else:
            self.labels = None
            self.labelgroup = None
            self.label_zorder = 15


    def view_prefs(self):
        """
        Create a table of the current map design elements (i.e. attributes)
            and their values.

        """

        pref_list = self.__dict__.keys()

        for pref, num in zip(pref_list, range(1, len(pref_list)+1)):
            print '\t', num,'. ', pref, ' : ', getattr(self, pref)

        print '\n'


    def edit_latlonlabels(self, lat_labels=None, lon_labels=None,
                          lat_labelspacing=None, lon_labelspacing=None,
                          latlon_fs=None):
        """
        """
        meridian_labels = [0,0,0,0]
        parallel_labels = [0,0,0,0]

        ind_dict = {'left' : 0,
                    'right' : 1,
                    'top' : 2,
                    'bottom' : 3}

        if lat_labels is not None:
            for la in lat_labels:
                parallel_labels[ind_dict[la]] = 1
            self.parallel_labels = parallel_labels
        if lon_labels is not None:
            for lo in lon_labels:
                meridian_labels[ind_dict[lo]] = 1
            self.meridian_labels = meridian_labels

        if lat_labelspacing is not None:
            self.latstep = lat_labelspacing

        if lon_labelspacing is not None:
            self.lonstep = lon_labelspacing

        if latlon_fs is not None:
            self.latlon_fs=latlon_fs



    def edit_resolution(self, resolution=None, area_threshold=None):
        """
        Adjust the map resolution and area threshold for plotting.
            Attributes will only be adjusted if not None.

        Keyword Arguments
        -----------------
        resolution : char
            Default 'c'.  ['c'|'l'|'i'|'h'|'f'].
            Crude, low, intermediate, high, full. The relative resolution of
            map boundaries.  Drops off by about 80 percent between datasets.
        area_threshold : int or string
            Default None.  ['auto'|1|10|100|1000|10000]
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
        Update the mapcorners.

        Parameters
        ----------
        mapcorners : list of floats
            Used to construct the map view for conic and cyl projections.
            Lower left longitude, latitude; upper right longitude, latitude.
        """

        self.mapcorners = mapcorners


    def edit_standard_pm(self, standard_pm):
        """
        Update the standard parallels and meridians.

        Parameters
        ----------
        standard_pm : list of floats
            For cyl and conic projections, the list creates standard parallels
                and meridians [lon_0, lat_0, lat_1, lat_2].
            For orthographic projection, lon_0 and lat_0 only are required.
                Sets the view from above the earth.
            For polar projections, lon_0 indicates the longitude that will be
                oriented N-S. lat_0 is replaced by the boundinglat,
                the lowest latitude that should appear on the map.
                lat_1 and lat_2 not required

        """

        self.standard_pm = standard_pm


    def swap_mapcolor(self):
        """
        Switch mapcolor from light to dark.

        """

        if self.mapcolor is 'light':
            self.mapcolor = 'dark'
        else:
            self.mapcolor = 'light'


    def edit_shapefiles(self, shp=None, delete_shp=None, overwrite_shps=False):
        """
        Edit the list of shapefiles.

        Keyword Arguments
        -----------------
        shp : tuple of list of tuples of strings
            Default None.  Tuples  are (File, color, linewidth)
            r'C:\programming\shapefiles\New_Shapefile'.
            Used to extend or overwrite current list of shapefiles.
        delete_shp : int or list of ints
            Default None.  The indices of shapefiles to delete.
        overwrite_shps : Boolean
            Default False.  Indicates whether to extend the current
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

        Keyword Arguments
        -----------------
        labelpath : string
            Default None.  Full or relative path to labelfile location
        labelgroup : string
            Default None.  ['all_1'|'all_2'|important'|'justcave'|'cave']
            The set of labels to apply.  'all_1' has caves instead of cities,
            'all_2' has cities instead of caves.

        """

        if labelpath is not None:
            self.labels = ml.labelfile_reader(labelpath)

        if labelgroup is not None:
            self.labelgroup = labelgroup

        if label_zorder is not None:
            self.label_zorder = label_zorder


    def clear_labels(self):
        """
        Reset label path and group information to None, label zorder to 15.

        """

        self.labels = None
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


    def make_basemap(self, figsize, ax=None):
        """
        Takes the MapDesign attributes plus a figure size and creates a map
            on which data can be plotted.

        Parameters
        ----------
        figsize : tuple of ints
            The size of the figure

        Keyword Arguments
        -----------------
        ax : axes instance
            Default None, figure and axis will be created.  Otherwise,
            basemap will be created on given axis.

        Returns
        -------
        fig : figure instance
            Returned only if ax=None and figure must be created.
        ax : axes instance
            The axis on which the basemap is drawn.
        cavemap : basemap instance
            A map ready for data plotting

        """

        if ax == None:
            # Create figure instance
            fig, ax = plt.subplots(1, 1, figsize=figsize)

        # Labels are left, right, top, bottom
        meridian_labels = self.meridian_labels
        parallel_labels = self.parallel_labels

        if self.projection is 'cconic':
            # Lambert conformal conic
            # Set area_threshold to 1000, to eliminate tiny lakes and islands
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

            meridian_labels = [0,0,0,0]
            parallel_labels = [0,0,0,0]

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
                                          self.labels, self.label_zorder)

        if len(self.shapefiles) > 0:
            # plot shapefiles
            for shp in self.shapefiles:
                cavemap.readshapefile(shp[0], 'hella_shapes', zorder=15,
                                      color=shp[1], linewidth=shp[2])

        # Draw countries, states, coastlines, and map boundary
        cavemap.drawcountries(zorder=self.zborder)
        cavemap.drawstates(zorder=14)
        cavemap.drawcoastlines(zorder=self.zborder, linewidth=1.5)

        # Draw and label meridians and parallels
        cavemap.drawmeridians(np.arange(-180, 180, self.long), labels=meridian_labels,
                                        fontsize=self.latlon_fs, zorder=11)
        cavemap.drawparallels(np.arange(-90, 90, 10), labels=parallel_labels,
                                        fontsize=self.latlon_fs, zorder=11)

        # Map color defaults to white for ortho projection
        if self.projection != 'ortho':
            if self.mapcolor == 'light':
                cavemap.drawmapboundary(zorder=16)
                cavemap.fillcontinents(color='0.95', zorder=12, alpha=0.85)
            else:
                cavemap.drawmapboundary(zorder=16, fill_color='0.3')
                cavemap.fillcontinents(color='0.5', lake_color='0.3',
                                       zorder=12, alpha=0.90)
        else:
            cavemap.fillcontinents(color='0.99', zorder=12, alpha=0.85)

        try:
            return fig, ax, cavemap
        except:
            return cavemap


def get_colormap(colormap):
    """
    Dictionary of available colormaps.

    Parameters
    ----------
    colormap : string
        Short name of colormap to retrieve
        ['jet'|'blues'|'anomaly'|'heat'|'earth']

    Returns
    -------
    colorscheme : colormap
        The colormap, ready for colormapping

    """

    colormap_dict = {'jet': plt.cm.jet, 'blues': plt.cm.Blues,
                     'anomaly': plt.cm.RdBu, 'heat': plt.cm.gist_heat_r,
                     'earth': plt.cm.gist_earth}

    colorscheme = colormap_dict[colormap]

    return colorscheme


def make_cbar(data, datamap, orientation='horizontal', cbar_size=(20, 1.0),
              divisions=5, reverse_cbar=False, cbar_label=None,
              tick_fs=16, label_fs=18):
    """
    Make a colorbar on the same axis as datamap.

    Parameters
    ----------
    data : matplotlib pyplot
        The plot for which a colorbar is needed
    datamap : Axis object
        The axis on which data is plotted

    Keyword Arguments
    -----------------
    orientation : string
        Default 'horizontal'.  ['horizontal'|'vertical'].  The location of the
        colorbar relative to the map.
    cbar_size : tuple of int, float
        Default (20, 1.0).  Colobar (aspect, shrink).  The H/W ratio of the
        colorbar, the fractional size of the colorbar.
    divisions : int
        Default 5.  The number of tick divisions on the colorbar
    reverse_cbar : Boolean
        Default False. If True, colorbar is flipped over short axis.
        Value-color mapping is unaffected.
    cbar_label : string
        Default None.  Colorbar label.
    tick_fs : int
        Default 16.  Font size of ticks
    label_fs : int
        Default 18.  Font size of cbar_label

    """

    # Initialize colorbar
    cbar = plt.colorbar(data, orientation=orientation, pad=.05,
                        aspect=cbar_size[0], shrink=cbar_size[1])

    # Make pretty
    edit_cbar(cbar, orientation, divisions, reverse_cbar, cbar_label, tick_fs,
              label_fs)


def make_cax_cbar(fig, rect, data, orientation='horizontal', divisions=5,
                  reverse_cbar=False, cbar_label=None, tick_fs=16, label_fs=18):
    """
    Make a colorbar on a new axis.

    Parameters
    ----------
    fig : figure instance
        The figure that needs a colorbar
    rect : list of floats
        The colorbar position and size.  [Distance from left, distance from
        bottom, size in x dimension, size in y dimension]

    Keyword Arguments
    -----------------
    orientation : string
        Default 'horizontal'.  ['horizontal'|'vertical'].  The orientation of
        the colormapping within in the colorbar.
    cbar_size : tuple of int, float
        Default (20, 1.0).  Colobar (aspect, shrink).  The H/W ratio of the
        colorbar, the fractional size of the colorbar.
    divisions : int
        Default 5.  The number of tick divisions on the colorbar.  If
        `divisions` is None, then the tick locator will be skipped.
    reverse_cbar : Boolean
        Default False. If True, colorbar is flipped over short axis.
        Value-color mapping is unaffected.
    cbar_label : string
        Default None.  Colorbar label.
    tick_fs : int
        Default 16.  Font size of ticks
    label_fs : int
        Default 18.  Font size of cbar_label

    Returns
    -------
    cax : matplotlib axes instance
        The new colorbar.  Remove using fig.delaxes(cax)

    """

    # Initialize cax and colorbar on cax
    cax = fig.add_axes(rect)
    cbar = fig.colorbar(data, cax=cax, orientation=orientation)

    # Make pretty
    edit_cbar(cbar, orientation, divisions, reverse_cbar, cbar_label, tick_fs,
              label_fs)

    return cax


def edit_cbar(cbar, orientation, divisions, reverse_cbar, cbar_label, tick_fs,
              label_fs):
    """
    Make the colorbar pretty.  Adjust fontsizes, add label, get a reasonable
        number of nice ticks, etc.

    Parameters
    ----------
    cbar : colorbar instance
        The colorbar created in make_cbar() or make_cax_cbar().
    orientation : string
        ['horizontal'|'vertical'].  The orientation of
        the colormapping within in the colorbar.
    cbar_size : tuple of int, float
        Colobar (aspect, shrink).  The H/W ratio of the
        colorbar, the fractional size of the colorbar.
    divisions : int
        The number of tick divisions on the colorbar
    reverse_cbar : Boolean
        If True, colorbar is flipped over short axis.
        Value-color mapping is unaffected.
    cbar_label : string
        Colorbar label.
    tick_fs : int
        Font size of ticks
    label_fs : int
        Font size of cbar_label

    """

    # Adjust ticks and tick labels
    if divisions is None:
        cbar.locator = tk.MaxNLocator(divisions, integer=False)
    cbar.ax.tick_params(labelsize=tick_fs)
    cbar.update_ticks()

    # Initialize dictionary
    rotation_dict = {'vertical' : (270, 24),
                     'horizontal' : (0, 10)}

    # Reverse colorbar
    if reverse_cbar:
        if orientation is 'horizontal':
            cbar.ax.invert_xaxis()
        else:
            cbar.ax.invert_yaxis()

    # Label colorbar
    if cbar_label is not None:
        rotation, labelpad = rotation_dict[orientation]
        cbar.set_label(cbar_label, labelpad=labelpad, fontsize=label_fs,
                       rotation=rotation)

    # Cbar will have lines through it if mappable's alpha < 1
    cbar.set_alpha(1)
    cbar.draw_all()