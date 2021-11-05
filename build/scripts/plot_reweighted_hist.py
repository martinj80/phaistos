# plot_reweighted_hist.py --- Plot reweighted 2d histogram over columns from an observable file
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

import sys
import muninn_weights

try:
    from plot_hist import *

    # Muninn includes
    from muninn.details.CanonicalAverager import CanonicalAverager
except:
    print >> sys.stderr, "Script files not found. Try running from build dir instead\n"
    raise


# Command line options
parser.add_option("--muninn-file", dest="muninn_log_filename", type="string",
                  help="Path to Muninn log file")
parser.add_option("--muninn-column", dest="muninn_column", default="0|sample_\d+_\d+_([-.\d]+)", type="string",
                  help="Column index (and optionally: |reg-exp-filter>lower-boundary<upper-boundary)")
parser.add_option("--beta", dest="beta", default=1.0, type="float",
                  help="Inverse temperature at which to reweight.")



if __name__ == "__main__":

    # Parse command line options
    (options, args) = parser.parse_args()

    if not options.muninn_log_filename:
        parser.error("Muninn log file not specified.")

    # Setup PDF output
    pp = PdfPages(options.output_file)

    # Create figure
    figsize = get_limits(options.figsize)
    fig = plt.figure(figsize=get_limits(options.figsize))

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
    # Get Muninn weights
    muninn_weights = muninn_weights.get_muninn_weights(options, data, input_filenames)

    # Add weights as column in data matrix
    # data = numpy.column_stack((data, muninn_weights))
    data = numpy.lib.recfunctions.append_fields(data, names='muninn', data=muninn_weights)

    # Enforce data ranges (if specified)
    data = restrict_to_range(options, data,
                              [(x_index, x_range),
                               (y_index, y_range)])
    
    # Renormalize weights column
    # weight_index = data.shape[1]-1
    # normalization_constant = sum(data[:,weight_index])
    # data[:,weight_index] /= normalization_constant
    normalization_constant = sum(data['muninn'])
    data['muninn'] /= normalization_constant

    # Calculate histogram
    bins = [options.bins, options.bins]
    if figsize != None:
        bins = [options.bins, options.bins*(figsize[1]/float(figsize[0]))]
    hist, x_edges, y_edges = calc_histogram(options, data, x_index, y_index, 
                                            bins=bins,
                                            weights=data['muninn'],
                                            range=[get_limits(options.x_lim),
                                                   get_limits(options.y_lim)])

    # Extract column names from first input file (if available)
    x_lab, y_lab = get_axis_names(options, input_filenames[0], 
                                  x_index, y_index)

    matplotlib.rc('font', **{'size':18})

    # Plot
    plot_histogram(options, fig, hist, x_lab, y_lab, x_edges, y_edges)

    # Close file
    pp.savefig(fig, bbox_inches='tight')
    pp.close()
