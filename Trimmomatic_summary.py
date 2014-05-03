__author__ = 'Kevin Rue-Albrecht'


# Module argparse allows to process the arguments from the command line
import argparse
# Module os allows to get the current working directory
import os
# Module sys allows to collect flags and arguments from the command line
import sys


def check_valid_reports(files_in):
    for file_in in files_in:
        if not os.path.isfile(file_in):
            print("Error: Report file was not found: %s" % file_in)
            sys.exit(2)


def get_input_reports(list_reports):
    files_in = []
    with open(list_reports) as filein:
        for line in filein:
            files_in.append(line.strip())
    # sanity check: make sure that each file listed in the input file exists
    check_valid_reports(files_in)
    return files_in


def get_report_statistics(input_report):
    statistics = []
    statistics_line = None
    with open(input_report) as filein:
        for line in filein:
            if line.startswith("Input Read Pairs"):
                statistics_line = line
    if not statistics_line:
        print('Could not find line starting with "Input Read Pairs" in file: %s' % input_report)
        sys.exit(2)
    #
    first_colon_index = statistics_line.index(':')
    both_index = statistics_line.index('Both')
    statistics.append(statistics_line[first_colon_index + 2:both_index - 1])
    #
    first_surviving_index = statistics_line.index('Surviving')
    first_bracket_index = statistics_line.index('(')
    statistics.append(statistics_line[first_surviving_index + 11:first_bracket_index - 1])
    #
    first_percent_index = statistics_line.index('%')
    statistics.append(statistics_line[first_bracket_index + 1:first_percent_index])
    #print("Test: surviving_both_percent: %s%%" % surviving_both_percent)
    second_surviving_index = statistics_line.index('Surviving', first_percent_index)
    second_bracket_index = statistics_line.index('(', second_surviving_index)
    statistics.append(statistics_line[second_surviving_index + 11:second_bracket_index - 1])
    #print("Test: surviving_forward_only_count: %s" % surviving_forward_only_count)
    second_percent_index = statistics_line.index('%', second_bracket_index)
    statistics.append(statistics_line[second_bracket_index + 1:second_percent_index])
    #print("Test: surviving_forward_only_percent: %s%%" % surviving_forward_only_percent)
    third_surviving_index = statistics_line.index('Surviving', second_percent_index)
    third_bracket_index = statistics_line.index('(', third_surviving_index)
    statistics.append(statistics_line[third_surviving_index + 11:third_bracket_index - 1])
    #print("Test: surviving_reverse_only_count: %s" % surviving_reverse_only_count)
    third_percent_index = statistics_line.index('%', third_bracket_index)
    statistics.append(statistics_line[third_bracket_index + 1:third_percent_index])
    #print("Test: surviving_reverse_only_percent: %s%%" % surviving_reverse_only_percent)
    dropped_index = statistics_line.index('Dropped', third_percent_index)
    fourth_bracket_index = statistics_line.index('(', dropped_index)
    statistics.append(statistics_line[dropped_index + 9:fourth_bracket_index - 1])
    #print("Test: dropped_count: %s" % dropped_count)
    fourth_percent_index = statistics_line.index('%', fourth_bracket_index)
    statistics.append(statistics_line[fourth_bracket_index + 1:fourth_percent_index])
    #print("Test: dropped_percent: %s%%" % dropped_percent)
    return statistics


def get_reports_statistics(input_reports):
    statistics = []
    for input_report in input_reports:
        statistics.append(get_report_statistics(input_report))
    return statistics


def write_report_statistics(report_statistics, output):
    with open(output, 'w') as fileout:
        fileout.write('\t'.join(['Input_count', 'Both_count', 'Both_percent', 'Forward_count', 'Forward_percent',
                                 'Reverse_count', 'Reverse_percent', 'Dropped_count', 'Dropped_percent\n']))
        for statistics in report_statistics:
            statistics.append('\n')
            fileout.write('\t'.join(statistics))


def __main__():
    # Informative message: print the current path
    print("Info: Working directory: {0:s}".format(os.getcwd()))
    # Define the argument parser
    parser = argparse.ArgumentParser(description="This package summarises multiple reports from Trimmomatic.")
    # Mandatory argument: the input filename
    parser.add_argument('-l', '--list-reports', required=True,
                        help="Name of the file listing the Trimmomatic report files. One file name per line",
                        metavar='reports.txt')
    # Mandatory argument: the output filename
    parser.add_argument('-o', '--output', required=True,
                        help="Name of the file to write the output to.",
                        metavar='summary.txt')
    # parse command line options according to the rules defined above
    args = parser.parse_args(sys.argv[1:])
    #print("Test: args: %s\n" % args)
    # sanity check: make sure that the input file exists
    if not os.path.isfile(args.list_reports):
        print("Error: File of reports was not found: %s" % args.list_reports)
        sys.exit(2)
    # Get the list of report files to summarise
    input_reports = get_input_reports(args.list_reports)
    #print("Test: input_reports: %s\n" % input_reports)
    # Collect the statistics from each individual report
    report_statistics = get_reports_statistics(input_reports)
    #print("Test: report_statistics: %s\n" % report_statistics)
    # Write the report statistics in the output file
    write_report_statistics(report_statistics, args.output)
    #
    print("Info: Successfully compiled trimmomatic results in file: %s" % args.output)


if __name__ == "__main__":
    __main__()

