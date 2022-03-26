# plot_columns.py --- Plot columns from an observable file
# Copyright (C) 2008-2013 Wouter Boomsma
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


import numpy
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages
import optparse
import re
import StringIO

# Command line options
usage = "usage: %prog [options] observable_file"
parser = optparse.OptionParser(usage=usage)
parser.add_option("--output", dest="output_file", default="plot.pdf",
                  help="Output filename")
parser.add_option("--figsize", dest="figsize", type="string",
                  help="Size of figure")
parser.add_option("--title", dest="title", default=None, type="string",
                  help="Plot title")
parser.add_option("--xlab", dest="x_lab", default=None, type="string",
                  help="label for x-axis")
parser.add_option("--ylab", dest="y_lab", default=None, type="string",
                  help="label for y-axis")
parser.add_option("--xlim", dest="x_lim", type="string",
                  help="Limits for x-axis data in format start:end")
parser.add_option("--ylim", dest="y_lim", type="string",
                  help="Limits for y-axis data in format start:end")


def get_legend_labels(options):
    'Setup legend labels'

    if options.legend_labels == None:
        return None
    
    return options.legend_labels.split(',')


def get_limits(option):
    'Setup limits for plot axis '

    limits = None
    if option:
        lower_str, upper_str = option.split(":")
        limits = [None,None]
        if lower_str:
            limits[0] = float(lower_str)
        if upper_str:
            limits[1] = float(upper_str)

    return limits


def setup_range(data, column_index, range):
    'Setup range (handles None values)'

    if range == None or (range[0]==None and range[1]==None):
        range = data[data.dtype.names[column_index]].min(), data[data.dtype.names[column_index]].max() 
        # range = data[:,column_index].min(), data[:,column_index].max() 
    elif range[0] == None:
        range = data[data.dtype.names[column_index]].min(), range[1] 
        # range = data[:,column_index].min(), range[1] 
    elif range[1] == None:
        range = range[0], data[data.dtype.names[column_index]].max()
        # range = range[0], data[:,column_index].max()
    return range


def get_column_info(index_str):
    'Extract index, filter and range options from a command line option string'

    # Create list of tokens
    prev_index = 0
    tokens = []
    for match in re.finditer("[|><]", index_str):
        tokens.append(index_str[prev_index:match.start()])
        prev_index = match.start()
    tokens.append(index_str[prev_index:])

    # Extract index, filter and range information
    index = None
    filter = ".*"
    range_lower = -float('inf')
    range_upper =  float('inf')
    for token in tokens:
        if token[0] == "|":
            filter = token[1:]
        elif token[0] == ">":
            range_lower = float(token[1:])
        elif token[0] == "<":
            range_upper = float(token[1:])
        else:
            index = int(token)

    return index, filter, [range_lower, range_upper]



def create_filter_converter(filter):
    'Create a converter function from a filter string'

    filter_re = re.compile(filter)
    converter = lambda x: float((filter_re.findall(x) + ['nan'])[0])
    return converter


def get_input_filenames(filenames):
    'Extract filenames from command line, also dealing with "-" special name'

    for i,filename in enumerate(filenames):
        if filename == '-':
            filenames[i] = StringIO.StringIO(sys.stdin.read())

    return filenames



def read_data(options, input_filename, 
              column_modifiers = []):
    'Extract data from input file'

    # Setup converters corresponding to filters
    converters = {}
    for column_modifier in column_modifiers:
        index, filter = column_modifier
        converters[index] = create_filter_converter(filter)

    filename = input_filename
    if isinstance(input_filename, StringIO.StringIO):
        filename = StringIO.StringIO(input_filename.getvalue())

    # Read in data
    data = numpy.genfromtxt(filename, dtype=None, 
                            converters=converters)
    return data


def restrict_to_range(options, data,
                       column_ranges = []):
    'Truncate data outside specified boundaries'

    # Filter with upper and lower limits
    for column_range in column_ranges:
        index, range = column_range
        # data = data[data[:,index]<range[1]]
        # data = data[data[:,index]>range[0]]
        data = data[data[data.dtype.names[index]]<range[1]]
        data = data[data[data.dtype.names[index]]>range[0]]

    return data


def get_axis_names(options, input_filename, x_index, y_index):
    'Attempt to extract column names from data file'

    if isinstance(input_filename, StringIO.StringIO):
        input_file = input_filename
    else:
        input_file = open(input_filename)
    header = input_file.readline()
    x_lab, y_lab = None,None
    if header[0] == '#':
        names = header[1:].strip().split()
        x_lab = names[x_index]
        y_lab = names[y_index]
    if options.x_lab:
        x_lab = options.x_lab
    if options.y_lab:
        y_lab = options.y_lab
    
    return x_lab, y_lab



if __name__ == "__main__":

    # Command line option
    parser.add_option("--using", dest="using", default="0|sample_(\d+)_:1", type="string",
                      help="Which columns to use separated by ':' (each field can optionally specify: |reg-exp-filter>lower-boundary<upper-boundary)")
    parser.add_option("--legend-labels", dest="legend_labels", default=None, type="string",
                      help="Labels to use when creating legend (ex: --legend-labels file1,file2)")
    (options, args) = parser.parse_args()

    # Read input filenames from arguments
    input_filenames = get_input_filenames(args)

    # Setup PDF output
    pp = PdfPages(options.output_file)

    # Create figure
    fig = plt.figure(figsize=get_limits(options.figsize))

    # Get filter indices and filters
    x_index_str, y_index_str = options.using.split(':')
    x_index, x_filter, x_range = get_column_info(x_index_str)
    y_index, y_filter, y_range = get_column_info(y_index_str)

    # Setup legends
    legend_labels = get_legend_labels(options)
    if legend_labels == None:
        legend_labels = input_filenames

    # Iterate over input filenames
    for i, input_filename in enumerate(input_filenames):

        # Extract data
        data = read_data(options, input_filename, 
                         [(x_index, x_filter),
                          (y_index, y_filter)])

        # Enforce data ranges (if specified)
        data = restrict_to_range(options, data,
                                 [(x_index, x_range),
                                  (y_index, y_range)])

        # Plot
        plt.plot(data[data.dtype.names[x_index]], data[data.dtype.names[y_index]],
                 label=legend_labels[i])


    # Extract column names from first input file (if available)
    x_lab, y_lab = get_axis_names(options, input_filenames[0], x_index, y_index)

    # Set axes labels 
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)

    # Set title
    if options.title:
        plt.title(options.title)

    # Set legend
    plt.legend(prop={'size':'small'})

    # Set axis limits
    plt.xlim(get_limits(options.x_lim))
    plt.ylim(get_limits(options.y_lim))

    # Close file
    pp.savefig(fig)
    pp.close()
