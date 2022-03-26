# muninn_weights.py --- Extract muninn weights from log file
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

import numpy.lib.recfunctions

try:
    from plot_hist import *

    # Muninn includes
    from muninn.details.CanonicalAverager import CanonicalAverager
except:
    print >> sys.stderr, "Script files not found. Try running from build dir instead\n"
    raise


def get_muninn_weights(options, data, input_filenames):
    'Calculate muninn weights'
    
    # Get Muninn column index and filter
    muninn_column_index, muninn_column_filter, muninn_column_limits = get_column_info(options.muninn_column)
    muninn_column_converter = create_filter_converter(muninn_column_filter)

    # Iterate over filenames and concatenate arrays containing muninn reaction coordinate
    muninn_column_data_list = []
    for input_filename in input_filenames:
        muninn_column_data = read_data(options, input_filename, 
                                       [(muninn_column_index, muninn_column_filter)])
        muninn_column_data_list.append(muninn_column_data)
    muninn_column_data = numpy.concatenate(muninn_column_data_list)

    # Instantiate Muninn Averager
    cp = CanonicalAverager(options.muninn_log_filename, -1)

    # Get weights
    # muninn_weights = numpy.array(cp.calc_weights(muninn_column_data[:,muninn_column_index], options.beta))
    muninn_weights = numpy.array(cp.calc_weights(muninn_column_data[muninn_column_data.dtype.names[muninn_column_index]], options.beta))

    return muninn_weights



if __name__ == "__main__":

    usage = "usage: %prog [options] observable_file"
    parser = optparse.OptionParser(usage=usage)

    # Command line options
    parser.add_option("--muninn-file", dest="muninn_log_filename", type="string", default="muninn.txt",
                      help="Path to Muninn log file")
    parser.add_option("--muninn-column", dest="muninn_column", default="0|sample_\d+_\d+_([-.\d]+)", type="string",
                      help="Column index (and optionally: |reg-exp-filter>lower-boundary<upper-boundary)")
    parser.add_option("--beta", dest="beta", default=1.0, type="float",
                      help="Inverse temperature at which to reweight.")
    parser.add_option("--output-weights-only", dest="output_weights_only", action="store_true", default=False, 
                      help="Whether only weight column should be outputted")

    # Parse command line options
    (options, args) = parser.parse_args()

    if not options.muninn_log_filename:
        parser.error("Muninn log file not specified.")

    # Read input filenames from arguments
    input_filenames = get_input_filenames(args)

    # Extract data
    data  = extract_merged_data(options, input_filenames)
    
    # Get Muninn weights
    muninn_weights = get_muninn_weights(options, data, input_filenames)

    # Add weights as column in data matrix
    # data = numpy.column_stack((data, muninn_weights))
    data = numpy.lib.recfunctions.append_fields(data, names='muninn', data=muninn_weights)

    # Enforce data ranges (if specified)
    data = restrict_to_range(options, data)

    # Renormalize weights column
    normalization_constant = sum(data['muninn'])
    data['muninn'] /= normalization_constant

    # print weights
    if options.output_weights_only:
        numpy.savetxt(sys.stdout, data['muninn'])
    else:
        numpy.savetxt(sys.stdout, data, delimiter=" ", fmt="%s\t"*len(data.dtype.names))

    
    
