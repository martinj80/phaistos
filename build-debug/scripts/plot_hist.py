# plot_hist.py --- Plot 2d histogram over columns from an observable file
# Copyright (C) 2008-2013 Wouter Boomsma, Jes Frellsen
#
# This file is part of Phaistos
#
# Phaistos is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Phaistos is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.


from plot_columns import *
import numpy.lib.recfunctions

# Command line options
parser.add_option("--zlim", dest="z_lim", type="string",
                  help="Limits for z-axis data in format start:end")
parser.add_option("--bins", dest="bins", default="100", type="int",
                  help="Number of bins to use for histogram")
parser.add_option("--using", dest="using", default="1:0|sample_\d+_\d+_([-.\d]+)", type="string",
                      help="Which columns to use separated by ':' (each field can optionally specify: |reg-exp-filter>lower-boundary<upper-boundary)")
parser.add_option("--interpolation", dest="interpolation", default="none", type="string",
                  help="Which type of interpolation to use (see help page for matplotlib.pyplot.imshow for details))")
parser.add_option("--log-mode", dest="log_mode", action="store_true", default=False, 
                  help="Plot log(P) values")
parser.add_option("--color-key", dest="color_key", action="store_true", default=False, 
                  help="Include color key")

def get_interpolation(options):
    'Extract interpolation type and radius from option value'

    interpolation_str = options.interpolation
    tokens = interpolation_str.split('(')

    interpolation = tokens[0]
    radius = 1
    if len(tokens) > 1:
       radius = float(tokens[1][:-1])
    return interpolation, radius


def extract_merged_data(options, input_filenames,
                        column_modifiers = []):
    'Extract data array constructed from multiple files'

    # Iterate over input filenames
    merged_data = []
    for input_filename in input_filenames:

        # Extract data
        data = read_data(options, input_filename, 
                         column_modifiers)

        # Add to list of arrays
        merged_data.append(data)

    # Concatenate arrays
    data = numpy.lib.recfunctions.stack_arrays(merged_data, autoconvert=True, usemask=False)

    return data


def calc_histogram(options, data, x_index, y_index,
                   bins, weights=None, range=None):
    'Calculate histogram data'

    if range != None:
        range[0] = setup_range(data, x_index, range[0])
        range[1] = setup_range(data, y_index, range[1])

    # hist, x_edges, y_edges = numpy.histogram2d(data[:,x_index], data[:,y_index], 
    hist, x_edges, y_edges = numpy.histogram2d(data[data.dtype.names[x_index]], data[data.dtype.names[y_index]], 
                                               bins=bins, 
                                               range=range,
                                               normed=True,
                                               weights=weights)
    if options.log_mode:
        # Take the logarithm (dealing with small values separately)
        limit = 1E-100
        iszero = hist<limit
        nonzero = hist>limit
        hist[iszero] = hist[nonzero].min()
        hist = numpy.log(hist)

    # Set histogram limits using z_lim option
    if options.z_lim:
        z_limits = get_limits(options.z_lim)
        if z_limits[0]:
            hist[hist<z_limits[0]] = z_limits[0]
        if z_limits[1]:
            hist[hist>z_limits[1]] = z_limits[1]        
    
    return hist, x_edges, y_edges


def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = numpy.mgrid[-size:size+1, -sizey:sizey+1]
    g = numpy.exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()


def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    from scipy import signal
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im,g, mode='same')
    return(improc)


def plot_histogram(options, pp, hist, 
                   x_lab, y_lab,
                   x_edges, y_edges):
    'Plot histogram'

    interpolation,interpolation_radius = get_interpolation(options)
    
    if interpolation == "blur":
        hist = blur_image(hist, 4)
        options.interpolation = 'none'

    # Set axis values
    extent = [x_edges[0], x_edges[-1], 
              y_edges[0], y_edges[-1]]

    # Plot
    plt.imshow(hist.T, extent=extent, aspect='auto', origin='lower', 
               interpolation=options.interpolation,
               filterrad=interpolation_radius)
        
    # Set axes labels 
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)

    # Set title
    if options.title:
        plt.title(options.title)

    # Create colorbar
    if options.color_key:
        plt.colorbar()


    


if __name__ == "__main__":

    # Parse command line options
    (options, args) = parser.parse_args()

    # Setup PDF output
    pp = PdfPages(options.output_file)

    # Create figure
    figsize = get_limits(options.figsize)
    fig = plt.figure(figsize=figsize)

    # Read input filenames from arguments
    input_filenames = get_input_filenames(args)

    # Get filter indices and filters
    x_index_str, y_index_str = options.using.split(':')
    x_index, x_filter, x_range = get_column_info(x_index_str)
    y_index, y_filter, y_range = get_column_info(y_index_str)

    # Extract data
    data = extract_merged_data(options, input_filenames,
                               [(x_index, x_filter),
                                (y_index, y_filter)])

    # Enforce data ranges (if specified)
    data = restrict_to_range(options, data,
                             [(x_index, x_range),
                              (y_index, y_range)])

    # Calculate histogram
    bins = [options.bins, options.bins]
    if figsize != None:
        bins = [options.bins, options.bins*(figsize[1]/float(figsize[0]))]
    hist, x_edges, y_edges = calc_histogram(options, data, x_index, y_index,
                                            bins=bins,
                                            range=[get_limits(options.x_lim),
                                                   get_limits(options.y_lim)])

    # Extract column names from first input file (if available)
    x_lab, y_lab = get_axis_names(options, input_filenames[0], 
                                  x_index, y_index)

    # Plot
    plot_histogram(options, fig, hist, x_lab, y_lab, x_edges, y_edges)
        
    # Close file
    pp.savefig(fig)
    pp.close()

