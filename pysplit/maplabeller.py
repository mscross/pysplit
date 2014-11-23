import os
from mpl_toolkits.basemap import Basemap


def map_labeller(cavemap, ax, labelgroups, label_list, labelzorder):
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
    labelgroups : string
        Determines what labelgroups are applied.
        ['all_1'|'all_2'|'cave'|'important'|'justcave']
    label_list : list of lists of strings or tuples of floats
        List of lists of labels and coordinates
        sea labels and coordinates
        cities
        oceans
        countries
        caves
    labelzorder : int
        Zorder of map labels.

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
                'zrd' : labelzorder,
                'off' : (0,0),
                'add' : ''}

    ocean_dict = {'name' : ocean_labels,
                  'coord' : ocean_coords,
                  'ha' : 'center',
                  'va' : 'center',
                  'fs' : 20,
                  'wt' : 'bold',
                  'fst' : 'italic',
                  'zrd' : labelzorder,
                  'off' : (0,0),
                  'add' : ''}

    country_dict = {'name' : country_labels,
                    'coord' : country_coords,
                    'ha' : 'center',
                    'va' : 'center',
                    'fs' : 20,
                    'wt' : 'bold',
                    'fst' : 'normal',
                    'zrd' : labelzorder,
                    'off' : (0,0),
                    'add' : ''}

    cave_dict = {'name' : cave_labels,
                 'coord' : cave_coords,
                 'ha' : 'right',
                 'va' : 'center',
                 'fs' : 15,
                 'wt' : 'normal',
                 'fst' : 'normal',
                 'zrd' : labelzorder,
                 'off' : (0.0, 1.0),
                 'add' : r'$\bigstar$'}

    city_dict = {'name' : cities_labels,
                 'coord' : cities_coords,
                 'ha' : 'right',
                 'va' : 'center',
                 'fs' : 15,
                 'wt' : 'normal',
                 'fst' : 'normal',
                 'zrd' : labelzorder,
                 'off' : (0.0, 1.0),
                 'add' : r'$\bullet$'}

    # Initialize dictionary of map label-group options
    label_masterdict = {'all_1' : [country_dict, ocean_dict, sea_dict, cave_dict],
                        'all_2' : [country_dict, ocean_dict, sea_dict, city_dict],
                        'important' : [country_dict, ocean_dict],
                        'city' : [city_dict, country_dict, ocean_dict],
                        'justcity': [city_dict],
                        'cave' : [cave_dict, country_dict, ocean_dict],
                        'justcave' : [cave_dict]}

    # Set option
    maplabels = label_masterdict[labelgroups]

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


def labelfile_reader(labelfile):
    """
    Opens and reads a text file with label information.

    If the file does not exist, the user will be prompted to:
    -Create a new label file at the given path, then continue attempting to
        label the map
    -Continue generating the map without labelgroups
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