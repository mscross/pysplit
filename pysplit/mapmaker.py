import numpy as np
import math
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
from mpl_toolkits.basemap import Basemap


def make_basemap(mapcorners, standard_pm, projection, figsize, labels,
                 shapefiles, mapcolor, ax=None):
    """
    Takes the mapview, projection type, further specifications as to labels and
        shapefiles and creates a basemap on which data can be plotted

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
    projection : string
        Indicates which projection to use.  Default is 'cyl'
            'cyl' : Equidistant cylindrical
            'eacyl'   : Equal Area cylindrical
            'cconic'  : Lambert Conformal Conic
            'eaconic' : Albers Equal Area Conic
            'ortho'   : Orthographic (globe)
            'npstere' : North polar steroegraphic (conformal)
            'spstere' : South polar steroegraphic (conformal)
            'nplaea'  : North polar azimuthal (equal area)
            'splaea'  : South polar azimuthal (equal area)
    figsize : tuple of ints
        The size of the figure
    labels : tuple of strings
        Default (None,None).  Label group, label file full or relative path.
        Determines what label groups are applied, if any:
        Label group choices:  ['all'|'important'|'justcave'|'cave']
        If path to label file does not exist, the user will be presented
        with several options, including one to make and edit a label file
    shapefiles : tuple of strings
        Default is (None, None, None).  File, color, linewidth
        r'C:\programming\shapefiles\New_Shapefile'
    mapcolor : string
        The map colorscheme (in grayscale). ['light'|'dark']
        'Dark' is not available for 'ortho' projections

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

    if projection is 'cconic':
        # Lambert conformal conic
        # Set area_threshold to 1000, to eliminate tiny lakes and islands
        cavemap = Basemap(llcrnrlon = mapcorners[0],
                          llcrnrlat = mapcorners[1],
                          urcrnrlon = mapcorners[2],
                          urcrnrlat = mapcorners[3],
                          projection = 'lcc',
                          lat_1 = standard_pm[2],
                          lat_2 = standard_pm[3],
                          lon_0 = standard_pm[0],
                          area_thresh = 1000,
                          resolution = 'h',
                          ax = ax)

    elif projection is 'eaconic':
        # Albers equal area conic
        cavemap = Basemap(llcrnrlon = mapcorners[0],
                          llcrnrlat = mapcorners[1],
                          urcrnrlon = mapcorners[2],
                          urcrnrlat = mapcorners[3],
                          projection = 'aea',
                          lat_1 = standard_pm[2],
                          lat_2 = standard_pm[3],
                          lon_0 = standard_pm[0],
                          lat_0 = standard_pm[1],
                          area_thresh = 1000,
                          resolution = 'h',
                          ax = ax)

    elif projection is 'eacyl':
        # equal area cylindrical
        cavemap = Basemap(llcrnrlon = mapcorners[0],
                          llcrnrlat = mapcorners[1],
                          urcrnrlon = mapcorners[2],
                          urcrnrlat = mapcorners[3],
                          projection = 'cea',
                          area_thresh = 1000,
                          resolution = 'h',
                          ax = ax)

    elif projection is 'ortho':
        # the earth
        cavemap = Basemap(projection = 'ortho',
                          lon_0 = standard_pm[0],
                          lat_0 = standard_pm[1],
                          area_thresh = 1000,
                          resolution = 'h',
                          ax = ax)

        meridian_labels = [0,0,0,0]
        parallel_labels = [0,0,0,0]

    elif projection[1:] == 'plaea' or projection[1:] == 'pstere':
        # Polar steroegraphic (conformal) or polar azimuthal(equal-area)
        cavemap = Basemap(projection = projection,
                          boundinglat = standard_pm[1],
                          lon_0 = standard_pm[0],
                          area_thresh = 1000,
                          resolution = 'h',
                          ax = ax)

    else:
        # Projection is 'cyl', Basemap's default equidist. cyl projection
        cavemap = Basemap(llcrnrlon = mapcorners[0],
                          llcrnrlat = mapcorners[1],
                          urcrnrlon = mapcorners[2],
                          urcrnrlat = mapcorners[3],
                          area_thresh = 1000,
                          resolution = 'h',
                          ax = ax)

    # Draw labels
    if labels[0] is not None:
        label_list = labelfile_reader(labels[1])
        if label_list is not None:
            cavemap = map_labeller(cavemap, ax, labels[0], label_list)

    if shapefiles[0] is not None:
        # plot shapefiles
        cavemap.readshapefile(shapefiles[0], 'hella_shapes', zorder=15,
                              color=shapefiles[1], linewidth=shapefiles[2])

    # Draw countries, states, coastlines, and map boundary
    cavemap.drawcountries(zorder=14)
    cavemap.drawstates(zorder=14)
    cavemap.drawcoastlines(zorder=14)

    # Draw and label meridians and parallels
    cavemap.drawmeridians(np.arange(-180, 180, 20), labels=meridian_labels,
                                    fontsize=20, zorder=11)
    cavemap.drawparallels(np.arange(-90, 90, 10), labels=parallel_labels,
                                    fontsize=20, zorder=11)

    if projection != 'ortho':
        if mapcolor == 'light':
            cavemap.drawmapboundary(zorder=16)
            cavemap.fillcontinents(color='0.95', zorder=12, alpha=0.90)
        else:
            cavemap.drawmapboundary(zorder=16, fill_color='0.3')
            cavemap.fillcontinents(color='0.5', lake_color='0.3',
                                   zorder=12, alpha=0.90)
    else:
        cavemap.fillcontinents(color='0.99', zorder=12, alpha=0.90)


    return cavemap



def labelfile_generator(labelfile):
    """
    Generates a label file template.

    Parameters
    ----------
    labelfile : string
        Full or relative path to template location

    """

    # Open new file
    labelfile = open(labelfile, 'w')

    labels = ['SEAS\n',
              '  16.00    88.50     Bay of\\nBengal\n',
              ' -15.05   115.00     South\\nChina\\nSea\n',
              '  27.00   125.00     East\\nChina\\nSea\n',
              'CITIES\n',
              '  39.91   116.39     Beijing\n',
              '  32.05   118.77     Nanjing\n',
              '  25.27   110.28     Guilin\n',
              'COUNTRIES\n',
              '  35.00   100.00     CHINA\n',
              'OCEANS\n',
              '  27.00   150.00     Pacific\\n\\nOcean\n',
              '  -5.00    70.00     Indian Ocean\n',
              'CAVES\n',
              '  32.29   119.05     Hulu\n',
              '  25.28   108.08     Dongge\n',
              '  39.41   115.39     Kulishu\n']

    labelfile.writelines(labels)
    labelfile.flush()
    labelfile.close()

    return None



def labelfile_reader(labelfile):
    """
    Opens and reads a text file with label information.

    If the file does not exist, the user will be prompted to:
    -Create a new label file at the given path, then continue attempting to
        label the map
    -Continue generating the map without labels
    -Enter a different path

    These three options will be presented until the user provides a valid
        labelfile or the user chooses to forgo labelling

    Parameters
    ----------
    labelfile : string
        Full or relative path to text file containing label information

    Returns
    -------
    labels : list of lists
        List of lists of labels and coordinates
        sea labels and coordinates
        cities
        oceans
        countries
        caves

    """

    labels=[]
    labeltypes = ['SEAS', 'CITIES', 'COUNTRIES', 'CAVES', 'OCEANS']

    # Look for file
    while True:
        if not os.path.exists(labelfile):
            print ('Label file does not exist.\n' +
                   'Press 1 to create a new label file\n'+
                   'Press 2 to make a map without labels\n'+
                   'Press 3 to enter a new label file path\n')
            useroption = raw_input()
            print type(useroption)
            if useroption == '1':
                labelfile_generator(labelfile)
                print 'Label file template generated at: '+ labelfile
                print 'Do not leave blank lines between labels/label headers'
                raw_input('Press enter when label file editing is done. ')
            elif useroption == '2':
                labels = None
                break
            elif useroption == '3':
                labelfile = raw_input('Enter new filename:\t')
        else:
            # Open file
            lfile = open(labelfile, 'r')

            # Read first line, a label header
            line = lfile.readline().strip()

            # Cycle through types of labels, breaking at end of file
            while True:
                if line =='':
                    break
                labellist=[]
                coordlist=[]

                # Run through labels within type, breaking at new types or EOF
                while True:
                    line = lfile.readline().strip()
                    if line in labeltypes or line == '':
                        break
                    latitude = float(line[:6])
                    longitude = float(line[6:14])
                    coord = (latitude, longitude)
                    place = line[19:]
                    place = place.replace('\\n', '\n')
                    labellist.append(place)
                    coordlist.append(coord)

                labels.append(labellist)
                labels.append(coordlist)

            break

    return labels



def map_labeller(cavemap, ax, labels, label_list):
    """
    Writes labels on a map.

    Use labelfile_generator() and labelfile_reader() to get label_list in the
        appropriate order and format.

    Parameters
    ----------
    cavemap : plot
        Instance of basemap to be labelled
    ax : axes instance
        axis of cavemap
    labels : string
        Determines what labels are applied.
        ['all'|'cave'|'important'|'justcave']
    label_list : list of lists of strings or tuples of floats
        List of lists of labels and coordinates
        sea labels and coordinates
        cities
        oceans
        countries
        caves

    Returns
    -------
    cavemap : plot
        Labelled instance of basemap

    """

    # Initialize lists of labels and coordinates
    sea_labels = label_list[0]
    sea_coords = label_list[1]

    cities_labels = label_list[2]
    cities_coords = label_list[3]

    ocean_labels = label_list[4]
    ocean_coords = label_list[5]

    country_labels = label_list[6]
    country_coords = label_list[7]

    cave_labels = label_list[8]
    cave_coords = label_list[9]


    # Initialize formatting dictionaries
    sea_dict = {'name' : sea_labels,
                'coord' : sea_coords,
                'ha' : 'center',
                'va' : 'center',
                'fs': 16,
                'wt' : 'normal',
                'fst' : 'italic',
                'zrd' : 15,
                'off' : (0,0),
                'add' : ''}

    ocean_dict = {'name' : ocean_labels,
                  'coord' : ocean_coords,
                  'ha' : 'center',
                  'va' : 'center',
                  'fs' : 20,
                  'wt' : 'bold',
                  'fst' : 'italic',
                  'zrd' : 15,
                  'off' : (0,0),
                  'add' : ''}

    country_dict = {'name' : country_labels,
                    'coord' : country_coords,
                    'ha' : 'center',
                    'va' : 'center',
                    'fs' : 20,
                    'wt' : 'bold',
                    'fst' : 'normal',
                    'zrd' : 15,
                    'off' : (0,0),
                    'add' : ''}

    cave_dict = {'name' : cave_labels,
                 'coord' : cave_coords,
                 'ha' : 'right',
                 'va' : 'center',
                 'fs' : 15,
                 'wt' : 'normal',
                 'fst' : 'normal',
                 'zrd' : 15,
                 'off' : (0.0, 1.0),
                 'add' : r'$\bigstar$'}

    # Initialize dictionary of map label-group options
    label_masterdict = {'all' : [country_dict, ocean_dict, sea_dict, cave_dict],
                        'important' : [country_dict, ocean_dict],
                        'cave' : [cave_dict, country_dict, ocean_dict],
                        'justcave' : [cave_dict]}

    # Set option
    maplabels = label_masterdict[labels]

    # Plot labels
    for j in maplabels:

        i = 0

        while i < len(j['name']):

            x, y = cavemap(j['coord'][i][1] + j['off'][1],
                           j['coord'][i][0] + j['off'][0])

            ax.text(x, y, j['name'][i] + j['add'],
                    horizontalalignment = j['ha'], verticalalignment = j['va'],
                    fontsize = j['fs'], weight = j['wt'], fontstyle = j['fst'],
                    zorder = j['zrd'])
            i = i + 1

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