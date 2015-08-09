from __future__ import division, print_function
import numpy as np


def grid_data(x, y, data, cell_value, binsize):
    """
    Place unevenly spaced 2D data on a grid by 2D binning using nearest
    neighbor interpolation.

    Originally Example 3 by ccampo (2010-07-11) from
    http://wiki.scipy.org/Cookbook/Matplotlib/Gridding_irregularly_spaced_data,

    Adapted for PySPLIT by mscross, 2014-06-12.

    Parameters
    ----------
    x : 1D ndarray of scalars
        The independent data x-axis of the grid (longitude)
    y : 1D ndarray of scalars
        The independent data y-axis of the grid (latitude)
    z : 1D ndarray of scalars
        The unevenly spaced dependent data.
        For example, specific humidity at each x,y point along a trajectory
    cell_value : string
        Determines the value of each cell from the contents of the bin.
        ['median'|'mean'|'cumulative'|'max'|'min'|'range'|'stdev']
    binsize : float
        The width, height of each bin.  Only square bins supported.

    Returns
    -------
    grid : masked 2D ndarray of scalars
        The evenly gridded data.  The value of each cell in relation to the
        bin contents is determined by ``cell_value``.
        Invalid values are masked.
    xi : 2D ndarray of floats
        The grid of x bin bounds
    yi : 2D ndarray of floats
        The grid of y bin bounds
    bins : 2D ndarray of floats
        A grid the same shape as ``grid``, except the value of each cell is
        the number of points in the bin.
    wherebin : 2D list
        A 2D list the same shape as ``grid`` and ``bins`` where each cell
        contains the indices of ``data`` that correspond to the values stored
        in the particular bin.

    """

    # Initialize dictionary
    cell_value_dict = {'cumulative': np.sum,
                       'mean': np.mean,
                       'median': np.median,
                       'max': np.max,
                       'min': np.min,
                       'stdev': np.std,
                       'range': np.ptp}

    # Get extreme longitudes and latitudes
    xmin = x.min()
    xmax = x.max()

    ymin = y.min()
    ymax = y.max()

    # Make coordinate arrays
    xi = np.arange(xmin, xmax + binsize, binsize)
    yi = np.arange(ymin, ymax + binsize, binsize)

    xi, yi = np.meshgrid(xi, yi)

    # Make `grid`
    grid = np.zeros(xi.shape, dtype=x.dtype)
    nrow, ncol = grid.shape

    # Set up `bins` and `wherebin` if either/both are called for
    bins = np.copy(grid)
    wherebin = np.copy(grid)
    wherebin = wherebin.tolist()

    # Fill in the grid
    for row in range(nrow):
        for col in range(ncol):

            # Get x, y coordinates of current position
            xc = xi[row][col]
            yc = yi[row][col]

            # Find the position(s) in the original data array that
            # xc, yc correspond to
            # Get absolute values of all items in x, y - xc, yc
            posx = np.abs(x - xc)
            posy = np.abs(y - yc)

            # Get boolean array of where condition is met in both x, y
            # ibin.size == x.size == y.size
            ibin = np.logical_and(posx < binsize / 2.0, posy < binsize / 2.0)

            # Get array of indices of data points that are in the current pos
            ind = np.where(ibin)[0]

            # Fill bin (True values in ibin will put data value into bin)
            bin = data[ibin]

            # Set cell values in grid
            if bin.size != 0:
                binval = cell_value_dict[cell_value](bin)
                grid[row, col] = binval
            else:
                grid[row, col] = -999.0

            # Update `wherebin` and `bins` if necessary
            wherebin[row][col] = ind
            bins[row, col] = bin.size

    # Mask 'invalid' entries.  PySPLIT follows the convention of
    # -999.0 as a fill value for invalid or missing data
    grid = np.ma.masked_less_equal(grid, -999.0)

    return grid, xi, yi, bins, wherebin
