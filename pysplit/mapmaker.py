import matplotlib.pyplot as plt
import matplotlib.ticker as tk


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

    return cbar


def make_cax_cbar(fig, rect, data, orientation='horizontal', divisions=5,
                  reverse_cbar=False, cbar_label=None, tick_fs=16,
                  label_fs=18):
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

    return cax, cbar


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
    if divisions is not None:
        cbar.locator = tk.MaxNLocator(divisions, integer=False)
    cbar.ax.tick_params(labelsize=tick_fs)
    cbar.update_ticks()

    # Initialize dictionary
    rotation_dict = {'vertical' : (270, 24, 'bottom'),
                     'horizontal' : (0, 10, 'baseline')}

    # Reverse colorbar
    if reverse_cbar:
        if orientation is 'horizontal':
            cbar.ax.invert_xaxis()
        else:
            cbar.ax.invert_yaxis()

    # Label colorbar
    if cbar_label is not None:
        rotation, labelpad, valign = rotation_dict[orientation]
        cbar.set_label(cbar_label, labelpad=labelpad, fontsize=label_fs,
                       rotation=rotation, verticalalignment=valign)

    # Cbar will have lines through it if mappable's alpha < 1
    cbar.set_alpha(1)
    cbar.draw_all()
