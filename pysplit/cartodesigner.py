from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

class CartoDesign(object):
    
    """Class for holding map design elements."""

    def __init__(self,xlabel_style, ylabel_style, mapcorners=(-95, 20, -25, -64), latlon_spacing=(10, 20),
                 draw_labels=True, ylabels_right = True, ylabels_left = True,
                 xlabels_top = True, xlabels_bottom = True,
                 background=True, bg_service='Stamen', bg_style='terrain', xsize=5, alpha_gl=0.3):
        """
        Initialize ``CartoDesign`` instance.

        Parameters
        ----------
		xlabel_style, ylabel_style = None
			A dictionary passed through to ax.text on x and y label creation for styling of the text labels.	
        mapcorners : list of floats
			(Lon_1, Lat_1, Lon_2, Lat_2)
            Lower left longitude, latitude; upper right longitude, latitude.
		latlonspacing : int or float
            Default (10, 20).  Degrees between plotted lines of latitude.
		draw_labels	: Boolean
			Default True.  Label gridlines like axis ticks, around the edge
		ylabels_right, ylabels_left, xlabels_top, xlabels_bottom : Boolean
			Default True.  Show labels around the edge.
		Background	: Boolean
			Default True.  Show background.
		bg_service 	: Services used by Cartopy:
			'Stamen':
				bg_style: 'terrain', 'toner', 'watercolor'
			'OSM'	:
				bg_style: Not used				
			'Google':
				bg_style: 'street', 'satellite'
			'Esri'	:
				bg_style: 'World_Topo_Map', 'World_Shaded_Relief', 'World_Imagery',
						  'NatGeo_World_Map', 'World_Physical_Map', 'World_Street_Map'
		xsize	: int or float
			Default 5. Zoom level for bg_service tiles.
		alpha_gl: int or float
			Default 0.3. Set transparency of gridlines.
        """
        self.mapcorners = mapcorners
        self.latlon_spacing = latlon_spacing
        self.background = background
        self.bg_service = bg_service
        self.bg_style = bg_style
        self.draw_labels = draw_labels
        self.ylabels_left = ylabels_left
        self.ylabels_right = ylabels_right
        self.xlabels_top = xlabels_top
        self.xlabels_bottom = xlabels_bottom
        self.xlabel_style = xlabel_style
        self.ylabel_style = ylabel_style
        self.xsize = xsize
        self.alpha_gl=alpha_gl

    def make_cartopy(self, projection, ax=None, figsize=(10, 10)):
        self.projection = projection
        if ax is None:
            # Create figure instance
            fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(projection=projection))
        Lat_1 = self.mapcorners[1]
        Lat_2 = self.mapcorners[3]
        Lon_1 = self.mapcorners[0]
        Lon_2 = self.mapcorners[2]
        xlocs_1 = ((round(Lon_1/10))-1)*10
        xlocs_2 = ((round(Lon_2/10))+1)*10
        ylocs_1 = ((round(Lat_1/10))-1)*10
        ylocs_2 = ((round(Lat_2/10))+1)*10		
        extent = [Lon_1, Lon_2, Lat_1, Lat_2]
        x_space = self.latlon_spacing[0]
        y_space = self.latlon_spacing[1]			
        xlocs=np.arange(xlocs_1, xlocs_2, x_space)
        ylocs=np.arange(ylocs_1, ylocs_2, y_space)
        import cartopy.io.img_tiles as cimgt
        if self.bg_service == 'Stamen':
            request = cimgt.Stamen(style=self.bg_style)
            ax.set_extent(extent, crs=self.projection)
            ax.add_image(request, self.xsize, interpolation='spline36')   		
        if self.bg_service == 'OSM':
            request = cimgt.OSM()
            ax.set_extent(extent, crs=self.projection)
            ax.add_image(request, self.xsize, interpolation='spline36')		
        if self.bg_service == 'Google':
            request = cimgt.GoogleTiles(style=self.bg_style)
            ax.set_extent(extent)
            ax.add_image(request, self.xsize, interpolation='spline36')
        if self.bg_service == 'Esri':
            request = cimgt.Esri(style=self.bg_style)
            ax.set_extent(extent, crs=self.projection)
            ax.add_image(request, self.xsize, interpolation='spline36')

        gl = ax.gridlines(draw_labels=self.draw_labels,xlocs=xlocs,ylocs=ylocs,alpha=self.alpha_gl)
        gl.ylabels_left = self.ylabels_left
        gl.ylabels_right = self.ylabels_right
        gl.xlabels_top = self.xlabels_top
        gl.xlabels_bottom = self.xlabels_bottom
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = self.xlabel_style
        gl.ylabel_style = self.ylabel_style  
        
        return ax