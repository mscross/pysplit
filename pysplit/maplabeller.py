from __future__ import division, print_function
import os


def map_labeller(basemap, which, labels, labelstyle, labelzorder):
    """
    Write labels on ``basemap``.

    Use ``labelfile_generator()`` and ``labelfile_reader()`` to get
    ``labels`` and ``labelstyle`` in the appropriate order and format.

    Parameters
    ----------
    basemap : ``Basemap`` instance
        Any ``Basemap`` instance.  For easy map creation, see ``MapDesign``
        class.
    which : list of strings
        Can include ['sea'|'country'|'ocean'|'place'|'city'].
        The labels to apply
    labels : list of lists of strings or tuples of floats
        List of lists of labels and coordinates
    labelzorder : int
        Zorder of map labels.
    labelstyle : list of dictionaries
        List of dictonaries that contain label formatting parameters

    Returns
    -------
    basemap : ``Basemap`` instance
        Labelled instance of ``Basemap``

    """
    # Initialize lists of labels and coordinates
    sea_labels = labels[0]
    sea_coords = labels[1]
    sea_style = labelstyle[0]

    cities_labels = labels[2]
    cities_coords = labels[3]
    cities_style = labelstyle[1]

    country_labels = labels[4]
    country_coords = labels[5]
    country_style = labelstyle[2]

    ocean_labels = labels[6]
    ocean_coords = labels[7]
    ocean_style = labelstyle[3]

    place_labels = labels[8]
    place_coords = labels[9]
    place_style = labelstyle[4]

    # Initialize formatting dictionaries
    sea_dict = {'place': sea_labels,
                'coord': sea_coords,
                'ha': 'center',
                'va': 'center',
                'fs': sea_style['fontsize'],
                'wt': sea_style['weight'],
                'fst': sea_style['fontstyle'],
                'zrd': labelzorder,
                'off': (0, 0),
                'add': ''}

    ocean_dict = {'place': ocean_labels,
                  'coord': ocean_coords,
                  'ha': 'center',
                  'va': 'center',
                  'fs': ocean_style['fontsize'],
                  'wt': ocean_style['weight'],
                  'fst': ocean_style['fontstyle'],
                  'zrd': labelzorder,
                  'off': (0, 0),
                  'add': ''}

    country_dict = {'place': country_labels,
                    'coord': country_coords,
                    'ha': 'center',
                    'va': 'center',
                    'fs': country_style['fontsize'],
                    'wt': country_style['weight'],
                    'fst': country_style['fontstyle'],
                    'zrd': labelzorder,
                    'off': (0, 0),
                    'add': ''}

    place_dict = {'place': place_labels,
                  'coord': place_coords,
                  'ha': 'right',
                  'va': 'center',
                  'fs': place_style['fontsize'],
                  'wt': place_style['weight'],
                  'fst': place_style['fontstyle'],
                  'zrd': labelzorder,
                  'off': (0.0, 1.0),
                  'add': r'$\bigstar$'}

    city_dict = {'place': cities_labels,
                 'coord': cities_coords,
                 'ha': 'right',
                 'va': 'center',
                 'fs': cities_style['fontsize'],
                 'wt': cities_style['weight'],
                 'fst': cities_style['fontstyle'],
                 'zrd': labelzorder,
                 'off': (0.0, 1.0),
                 'add': r'$\bullet$'}

    label_dict = {'sea': sea_dict,
                  'city': city_dict,
                  'place': place_dict,
                  'ocean': ocean_dict,
                  'country': country_dict}

    maplabels = [label_dict[label] for label in which]

    # Plot labels
    for j in maplabels:

        for i in range(0, len(j['place'])):

            x, y = basemap(j['coord'][i][1] + j['off'][1],
                           j['coord'][i][0] + j['off'][0])

            basemap.ax.text(x, y, j['place'][i] + j['add'],
                            horizontalalignment=j['ha'],
                            verticalalignment=j['va'],
                            fontsize=j['fs'], weight=j['wt'],
                            fontstyle=j['fst'], zorder=j['zrd'])

    return basemap


def labelfile_reader(labelfile):
    """
    Open and read a text file with label information.

    Parameters
    ----------
    labelfile : string
        Full or relative path to text file containing label information

    Returns
    -------
    labels : list of lists
        List of lists of labels and coordinates
    labelstyle : list of dictionaries
        List of dictonaries that contain label formatting parameters

    """
    labels = []
    labelstyle = []
    labeltypes = ['SEA', 'CITY', 'COUNTRY', 'PLACE', 'OCEAN']

    # Look for file
    if not os.path.exists(labelfile):
        raise OSError(labelfile, ' does not exist.\n',
                      'Generate a template at this location',
                      'using labelfile_generator(', labelfile, ')')

    # Open file
    with open(labelfile, 'r') as lfile:

        # Read first line, a label header
        line = lfile.readline().strip()

        # Cycle through types of labels, breaking at end of file
        while True:
            if line == '':
                break
            labellist = []
            coordlist = []

            # Run through labels within type
            while True:
                line = lfile.readline().strip()

                # If line is a new type
                if line in labeltypes:
                    # Get the next line, the style parameters
                    line = lfile.readline().strip()

                    if line == '':
                        break

                    # Break into substrings, one parameter per string
                    line = line.split('  ', 2)

                    # Split again over the '=', arrange into dictionary
                    for i in range(0, len(line)):
                        line[i] = line[i].strip().split('=')

                    styledict = {k: v for k, v in line}
                    labelstyle.append(styledict)

                    break

                elif line == '':
                    break

                else:
                    latitude = float(line[:6])
                    longitude = float(line[6:14])

                    coord = (latitude, longitude)

                    place = line[19:]
                    place = place.replace('\\n', '\n')

                    labellist.append(place)
                    coordlist.append(coord)

            if len(labellist) > 0:
                labels.append(labellist)
                labels.append(coordlist)

    return labels, labelstyle


def labelfile_generator(labelfile):
    """
    Generate a label file template.

    Parameters
    ----------
    labelfile : string
        Full or relative path to template location

    """
    # Open new file
    with open(labelfile, 'w') as labelfile:

        labels = ['File header\n',
                  'SEA\n',
                  '  fontstyle=italic   weight=normal   fontsize=16\n',
                  '  16.00    88.50     Bay of\\nBengal\n',
                  ' -15.05   115.00     South\\nChina\\nSea\n',
                  '  27.00   125.00     East\\nChina\\nSea\n',
                  'CITY\n',
                  '  fontstyle=normal   weight=normal   fontsize=15\n',
                  '  39.91   116.39     Beijing\n',
                  '  32.05   118.77     Nanjing\n',
                  '  25.27   110.28     Guilin\n',
                  'COUNTRY\n',
                  '  fontstyle=normal   weight=bold     fontsize=20\n',
                  '  35.00   100.00     CHINA\n',
                  'OCEAN\n',
                  '  fontstyle=italic   weight=bold     fontsize=20\n',
                  '  27.00   150.00     Pacific\\n\\nOcean\n',
                  '  -5.00    70.00     Indian Ocean\n',
                  'PLACE\n',
                  '  fontstyle=normal   weight=normal   fontsize=15\n',
                  '  32.29   119.05     Hulu\n',
                  '  25.28   108.08     Dongge\n',
                  '  39.41   115.39     Kulishu\n']

        labelfile.writelines(labels)
        labelfile.flush()
