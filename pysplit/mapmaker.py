import matplotlib.pyplot as plt
import matplotlib.ticker as tk


def make_cbar(data, ax, orientation='horizontal', cbar_size=(20, 1.0),
              reverse_cbar=False, **kwargs):
    """
    Make a colorbar on the same axis as ax.

    Parameters
    ----------
    data : PathCollection
        The mappable
    ax : Axis object
        The axis on which data is plotted

    Keyword Arguments
    -----------------
    orientation : string
        Default 'horizontal'.  ['horizontal'|'vertical'].  The location of the
        colorbar relative to the map.
    cbar_size : tuple of int, float
        Default (20, 1.0).  Colobar (aspect, shrink).  The H/W ratio of the
        colorbar, the fractional size of the colorbar.
    reverse_cbar : Boolean
        Default False. If True, colorbar is flipped over short axis.
        Value-color mapping is unaffected.

    Other Parameters
    ----------------
    kwargs : passed to edit_cbar()

    Returns
    -------
    cbar : matplotlib colorbar instance
        The new colorbar

    """

    # Initialize colorbar
    cbar = plt.colorbar(data, ax=ax, orientation=orientation,
                        aspect=cbar_size[0], shrink=cbar_size[1])

    # Reverse colorbar
    if reverse_cbar:
        if orientation is 'horizontal':
            cbar.ax.invert_xaxis()
        else:
            cbar.ax.invert_yaxis()

    # Make pretty
    edit_cbar(cbar, **kwargs)

    return cbar


def make_cax_cbar(fig, rect, data, orientation='horizontal',
                  reverse_cbar=False, extend='neither', **kwargs):
    """
    Make a colorbar on a new axis.

    Parameters
    ----------
    fig : figure instance
        The figure that needs a colorbar
    rect : list of floats
        The colorbar position and size.  [Distance from left, distance from
        bottom, size in x dimension, size in y dimension]
    data : PathCollection
        Mappable

    Keyword Arguments
    -----------------
    orientation : string
        Default 'horizontal'.  ['horizontal'|'vertical'].  The orientation of
        the colormapping within in the colorbar.
    cbar_size : tuple of int, float
        Default (20, 1.0).  Colobar (aspect, shrink).  The H/W ratio of the
        colorbar, the fractional size of the colorbar.
    reverse_cbar : Boolean
        Default False. If True, colorbar is flipped over short axis.
        Value-color mapping is unaffected.
    extend : string
        Default 'neither'.  ['both'|'neither'|'under'|'over'].
        Extend colorbar with pointed ends.

    Other Parameters
    ----------------
    kwargs : passed to edit_cbar()

    Returns
    -------
    cax : matplotlib axes instance
        The axis of the new colorbar.  Remove using fig.delaxes(cax)
    cbar : matplotlib colorbar instance
        The new colorbar

    """

    # Initialize cax and colorbar on cax
    cax = fig.add_axes(rect)
    cbar = fig.colorbar(data, cax=cax, orientation=orientation,
                        extend=extend)
    # Reverse colorbar
    if reverse_cbar:
        if orientation is 'horizontal':
            cbar.ax.invert_xaxis()
        else:
            cbar.ax.invert_yaxis()

    # Make pretty
    edit_cbar(cbar, **kwargs)

    return cax, cbar


def edit_cbar(cbar, divisions=5, cbar_label=None, tick_fs=16, label_fs=18,
              labelpad=24, rotation=0, tick_dir='out', tick_dim=(4, 2)):
    """
    Make the colorbar pretty.  Adjust fontsizes, add label, get a reasonable
        number of nice ticks, etc.

    Parameters
    ----------
    cbar : colorbar instance
        The colorbar created in make_cbar() or make_cax_cbar().

    Keyword Arguments
    -----------------
    divisions : int
        Default 5.  The number of nice ticks on the colorbar.  May be None.
    cbar_label : string
        Default None.  Colorbar label.
    tick_fs : int
        Default 16.  Font size of ticks
    label_fs : int
        Default 18.  Font size of cbar_label
    labelpad : int
        Default 24.  Spacing between tick labels and colorbar label
    rotation : int
        Default 0.  Rotation in degrees of label.
    tick_dir : string
        Default 'out'.  ['out'|'in'|'inout']
        Direction that ticks are pointing relative to colorbar
    tick_dim : tuple of floats
        Default (4, 2).  The (length, width) of the colorbar ticks

    """

    # Adjust ticks and tick labels
    if divisions is not None:
        cbar.locator = tk.MaxNLocator(divisions, integer=False)

    cbar.ax.tick_params(labelsize=tick_fs, direction=tick_dir,
                        length=tick_dim[0], width=tick_dim[1])
    cbar.update_ticks()

    # Label colorbar
    if cbar_label is not None:
        cbar.set_label(cbar_label, labelpad=labelpad, fontsize=label_fs,
                       rotation=rotation)

    # Cbar will have lines through it if mappable's alpha < 1
    cbar.set_alpha(1)
    cbar.draw_all()
