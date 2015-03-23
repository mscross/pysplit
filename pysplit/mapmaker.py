import matplotlib.pyplot as plt
import matplotlib.ticker as tk
import matplotlib.colors as clr
import numpy as np


def traj_scatter(data, lons, lats, cavemap, zorder=19, colormap=plt.cm.Blues,
                 edgecolor='none', size=25, sizedata=None, cnormalize=None,
                 snormalize=None, vmin=None, vmax=None, steps=11, **kwargs):
    """
    Scatter-plot of trajectory, trajectory group/cluster data.

    Parameters
    ----------
    data : 1D ndarray of floats, ints
        The data to plot as color change.
    lons : 1D ndarray of floats, ints
        `data` longitudes
    lats : 1D ndarray of floats, ints
        `data` latitudes
    cavemap : Basemap instance
        Initialize a basemap first using MapDesign.make_basemap()

    Keyword Arguments
    -----------------
    zorder : int
        Default 19.  Data zorder.
    colormap : colormap
        Default `plt.cm.Blues`.  Any matplotlib colormap.
    edgecolor : string, tuple
        Default `none`.  Any matplotlib-accepted color
    size : int
        Default 25.  Point size of data unless `sizedata` specified.  Then,
        will be multiplied with `sizedata`.
    sizedata : 1D ndarray of floats
        Default None.  The data to plot as a change in marker size.
    cnormalize : string
        Default None.  [None|'boundary'|'log'|'ln'|'sqrt']
        Normalization of color scale.  If 'boundary', will create a discrete
        color map with `steps` number of colors.  For other norms, colorbar
        ticklabels will be updated with 'log' but not 'ln', 'sqrt', because
        no corresponding matplotlib colors Normalize classes are available.
    snormalize : string
        Default None.  [None|'log'|'ln'|'sqrt'].  Similar to cnormalize,
        except 'boundary' not available and does not use Normalize.
    vmin : int or float
        Default None.
    vmax : int or float
    steps : int
        Only used in BoundaryNorm

    Other Parameters
    ----------------
    kwargs : passed to ax.scatter()

    Returns
    -------

    """

    cnormalize = str.lower(cnormalize)
    norm = None
    msg = ('Use `cbar.ax.set_yticklabels()` ' +
           'or cbar.ax.set_xticklabels()` to change tick labels')

    transform_dict = {'sqrt' : np.sqrt,
                      'log'  : np.log10,
                      'ln'   : np.ln}

    if cnormalize is 'boundary':
        if vmin is None:
            vmin = data.min()
        if vmax is None:
            vmax = data.max()
        bounds = np.linspace(vmin, vmax, steps)
        norm = clr.BoundaryNorm(bounds, colormap.N)
    elif cnormalize is 'log':
        norm = clr.LogNorm()
    elif cnormalize is 'ln':
        data = np.log(data)
        print msg, '\nnatural log normalization'
    elif cnormalize is 'sqrt':
        data = np.sqrt(data)
        print msg, '\nsqrt normalization'

    if sizedata is not None:
        if snormalize is not None:
            sizedata = transform_dict[snormalize](sizedata)
        size = sizedata * size

    cm = cavemap.scatter(lons, lats, c=data, s=size, cmap=colormap,
                         vmin=vmin, vmax=vmax, zorder=zorder,
                         edgecolor=edgecolor, norm=norm, latlon=True, **kwargs)

    return cavemap, cm


def traj_path(cavemap, lons, lats, color, lw, marker=None, linestyle='-',
              markeredgecolor='none', zorder=19, **kwargs):

    """
    Line plot of trajectory or cluster path

    Parameters
    ----------

    Keyword Arguments
    -----------------

    Other Parameters
    ----------------

    Returns
    -------

    """

    cavemap.plot(lons, lats, color, linewidth=lw, linestyle=linestyle,
                 marker=marker, latlon=True, zorder=zorder,
                 markeredgecolor=markeredgecolor, **kwargs)

    return cavemap


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


def random_colors(number_ofcolors):
    """
    """

    color_tmp = np.random.rand(number_ofcolors, 3)
    color_tmp = np.vsplit(color_tmp, number_ofcolors)
    colors = []
    for c in color_tmp:
        colors.append(c[0])
