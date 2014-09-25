import numpy as np
import math
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
from mpl_toolkits.basemap import Basemap
import maplabeller as ml


class MapDesign(object):
    """
    Class for holding map design elements

    """

    def __init__(self, mapcorners, standard_pm, projection='cyl',
                 mapcolor='light', shapefiles=[],
                 labels=None):
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
                'eacyl'   : Equal Area cylindrical
                'cconic'  : Lambert Conformal Conic
                'eaconic' : Albers Equal Area Conic
                'ortho'   : Orthographic (globe)
                'npstere' : North polar steroegraphic (conformal)
                'spstere' : South polar steroegraphic (conformal)
                'nplaea'  : North polar azimuthal (equal area)
                'splaea'  : South polar azimuthal (equal area)
        mapcolor : string
            Default 'light' The map grayscheme. ['light'|'dark']
            'Dark' is not available for 'ortho' projections
        shapefiles : list of tuples of strings
            Default is [].  (File, color, linewidth)
            r'C:\programming\shapefiles\New_Shapefile'
        labels : tuple of strings
            Default None.  Label group, label file full or relative path.
            Determines what label groups are applied, if any:
            Label group choices:  ['all'|'important'|'justcave'|'cave']
            If path to label file does not exist, the user will be presented
            with several options, including one to make and edit a label file

        """

        self.mapcorners = mapcorners
        self.standard_pm = standard_pm
        self.mapcolor = mapcolor
        self.shapefiles = shapefiles

        self.edit_projection(projection)

        if labels is not None:
            self.labels = ml.labelfile_reader(labels[1])
            self.labelgroup = labels[0]
        else:
            self.labels = None
            self.labelgroup = None


    def __dir__(self):
        """
        Get a list of the MapDesign attribute names

        """

        return [self.mapcorners, self.projection, self.standard_pm,
                self.mapcolor, self.shapefiles, self.labels, self.labelgroup]


    def view_prefs(self):
        """
        Prints the current design elements

        """
        pref_list = self.__dir__()

        for pref, num in zip(pref_list, range(1, len(pref_list)+1)):
            print '\t', num, '. ', pref, ' : ', getattr(self, pref)

        print '\n'


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
        Switch mapcolor from light to dark

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
        shp
        delete_shp
        overwrite_shps

        """

        if delete_shp is not None:
            del(self.shapefiles[delete_shp])

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
        projection

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
            print 'Projection not recognized.'
            print 'Choose a projection from the left column:'
            for proj in available_proj:
                print '\t', proj, '\t', available_proj[proj]


    def edit_labels(self, labels):

        if labels is not None:
            self.labels = ml.labelfile_reader(labels[1])
            self.labelgroup = labels[0]
        else:
            self.labels = None
            self.labelgroup = None


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
            Default None, figure will be created

        Returns
        -------
        cavemap : plot
           Basemap instance on which data can be plotted

        """

        if ax == None:
            # Create figure instance
            fig, ax = plt.subplots(1, 1, figsize=figsize)

        # Labels are left, right, top, bottom
        meridian_labels = [0,0,1,0]
        parallel_labels = [1,0,0,0]

        if self.projection is 'cconic':
            # Lambert conformal conic
            # Set area_threshold to 1000, to eliminate tiny lakes and islands
            cavemap = Basemap(llcrnrlon = self.mapcorners[0],
                              llcrnrlat = self.mapcorners[1],
                              urcrnrlon = self.mapcorners[2],
                              urcrnrlat = self.mapcorners[3],
                              projection = 'lcc',
                              lat_1 = self.standard_pm[2],
                              lat_2 = self.standard_pm[3],
                              lon_0 = self.standard_pm[0],
                              area_thresh = 1000,
                              resolution = 'h',
                              ax = ax)

        elif self.projection is 'eaconic':
            # Albers equal area conic
            cavemap = Basemap(llcrnrlon = self.mapcorners[0],
                              llcrnrlat = self.mapcorners[1],
                              urcrnrlon = self.mapcorners[2],
                              urcrnrlat = self.mapcorners[3],
                              projection = 'aea',
                              lat_1 = self.standard_pm[2],
                              lat_2 = self.standard_pm[3],
                              lon_0 = self.standard_pm[0],
                              lat_0 = self.standard_pm[1],
                              area_thresh = 1000,
                              resolution = 'h',
                              ax = ax)

        elif self.projection is 'eacyl':
            # equal area cylindrical
            cavemap = Basemap(llcrnrlon = self.mapcorners[0],
                              llcrnrlat = self.mapcorners[1],
                              urcrnrlon = self.mapcorners[2],
                              urcrnrlat = self.mapcorners[3],
                              projection = 'cea',
                              area_thresh = 1000,
                              resolution = 'h',
                              ax = ax)

        elif self.projection is 'ortho':
            # the earth
            cavemap = Basemap(projection = 'ortho',
                              lon_0 = self.standard_pm[0],
                              lat_0 = self.standard_pm[1],
                              area_thresh = 1000,
                              resolution = 'h',
                              ax = ax)

            meridian_labels = [0,0,0,0]
            parallel_labels = [0,0,0,0]

        elif projection[1:] == 'plaea' or projection[1:] == 'pstere':
            # Polar steroegraphic (conformal) or polar azimuthal(equal-area)
            cavemap = Basemap(projection = self.projection,
                              boundinglat = self.standard_pm[1],
                              lon_0 = self.standard_pm[0],
                              area_thresh = 1000,
                              resolution = 'h',
                              ax = ax)

        else:
            # Projection is 'cyl', Basemap's default equidist. cyl projection
            cavemap = Basemap(llcrnrlon = self.mapcorners[0],
                              llcrnrlat = self.mapcorners[1],
                              urcrnrlon = self.mapcorners[2],
                              urcrnrlat = self.mapcorners[3],
                              area_thresh = 1000,
                              resolution = 'h',
                              ax = ax)

        # Draw labels
        if self.labels is not None:
                cavemap = ml.map_labeller(cavemap, ax, self.labelgroup,
                                          self.labels)

        if len(self.shapefiles) > 0:
            # plot shapefiles
            for shp in self.shapefiles:
                cavemap.readshapefile(shp[0], 'hella_shapes', zorder=15,
                                      color=shp[1], linewidth=shp[2])

        # Draw countries, states, coastlines, and map boundary
        cavemap.drawcountries(zorder=14)
        cavemap.drawstates(zorder=14)
        cavemap.drawcoastlines(zorder=14)

        # Draw and label meridians and parallels
        cavemap.drawmeridians(np.arange(-180, 180, 20), labels=meridian_labels,
                                        fontsize=20, zorder=11)
        cavemap.drawparallels(np.arange(-90, 90, 10), labels=parallel_labels,
                                        fontsize=20, zorder=11)

        if self.projection != 'ortho':
            if self.mapcolor == 'light':
                cavemap.drawmapboundary(zorder=16)
                cavemap.fillcontinents(color='0.95', zorder=12, alpha=0.90)
            else:
                cavemap.drawmapboundary(zorder=16, fill_color='0.3')
                cavemap.fillcontinents(color='0.5', lake_color='0.3',
                                       zorder=12, alpha=0.90)
        else:
            cavemap.fillcontinents(color='0.99', zorder=12, alpha=0.90)


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