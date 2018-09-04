#/usr/bin/env python

from __future__ import division

__author__ = "Jesse Zaneveld"
__copyright__ = "Copyright 2016, The Karenina Project"
__credits__ = ["Jesse Zaneveld"]
__license__ = "GPL"
__version__ = "0.0.1-dev"
__maintainer__ = "Jesse Zaneveld"
__email__ = "zaneveld@gmail.com"
__status__ = "Development"

import karenina.visualization
from karenina.experiment import Experiment
from optparse import OptionParser
from optparse import OptionGroup
from os.path import join,isdir,realpath,dirname
import os
from pkg_resources import resource_filename
from os import makedirs
import pandas as pd


def make_option_parser():
    """Return an optparse OptionParser object"""

    parser = OptionParser(usage = "%prog -o ./simulation_results",
    description = "This script simulates microbiome " +
    "change over time using Ornstein-Uhlenbeck (OU) models.  These are " +
    "similar to Brownian motion models, with the exception that they " +
    "include reversion to a mean. Output is a tab-delimited data table " +
    "and figures.",
    version = __version__)

    required_options = OptionGroup(parser, "Required options")

    required_options.add_option('-o','--output', type="string",
    help='the output folder for the simulation results')

    parser.add_option_group(required_options)

    optional_options = OptionGroup(parser, "Optional options")

    optional_options.add_option('--pert_file_path',\
      default = os.path.abspath(resource_filename('karenina.data','set_xyz_lambda_zero.tsv')),\
      type = "string",\
      help = 'file path to a perturbation file specifying parameters for' +
      ' the simulation results [default: %default]')

    optional_options.add_option('--treatment_names',\
    default="control,destabilizing_treatment",type="string",\
    help="Comma seperated list of treatment named [default:%default]")

    optional_options.add_option('-n','--n_individuals',\
    default="35,35",type="string",\
    help='Comma-separated number of individuals to simulate per treatment.'+\
    'Note: This value must be enclosed in quotes. Example: "35,35". [default: %default]')

    optional_options.add_option('-t', '--n_timepoints',default=10, type="int",\
    help='Number of timepoints to simulate. (One number, which is the ' +
    'same for all treatments) [default: %default]')

    optional_options.add_option('-p','--perturbation_timepoint',\
    default=5,type="int",\
    help='Timepoint at which to apply a perturbation. Must be less than ' +
    '--n_timepoints [default: %default]')

    optional_options.add_option('-d','--perturbation_duration',\
    default=100,type="int",\
    help='Duration that the perturbation lasts. [default: %default]')

    optional_options.add_option('--interindividual_variation',
    default=0.01,type="float",help='Starting variability between ' +
    'individuals. [default: %default]')

    optional_options.add_option('--delta',default=0.25,type="float",
    help='Starting delta parameter for Brownian motion and ' +
    'Ornstein-Uhlenbeck processes. A higher number indicates more ' +
    'variability over time. [default: %default]')

    optional_options.add_option('-l','--L',default=0.20,type="float",
    help='Starting lambda parameter for Ornstein-Uhlenbeck processes. A ' +
    'higher number indicates a greater tendancy to revert to the mean ' +
    'value. [default: %default]')

    optional_options.add_option('--fixed_start_pos',default=None,type="string",
    help='Starting x,y,z position for all points, as comma separated ' +
    'floating point values, e.g. 0.0,0.1,0.2. If not supplied, starting ' +
    'positions will be randomized based on the interindividual_variation ' +
    'parameter [default: %default]')

    optional_options.add_option('-v','--verbose', action="store_true", dest="verbose", default=False,
                                help='-v, allows for verbose output' +
                                     ' [default: %default]')

    parser.add_option_group(optional_options)

    return parser


def check_perturbation_timepoint(perturbation_timepoint,n_timepoints):
    """
    Raise ValueError if perturbation_timepoint is < 0 or >n_timepoints

    :param perturbation_timepoint: defined timepoint for perturbation application
    :param n_timepoints: number of timepoints
    """

    if perturbation_timepoint and perturbation_timepoint >= n_timepoints:
        raise ValueError("Perturbation timepoint must be before the last timepoint")
    if perturbation_timepoint < 0:
        raise ValueError("Perturbation timepoint must be positive")


def ensure_exists(output_dir):
    """
    Ensure that output_dir exists

    :param output_dir: path to output directory
    """
    try:
        makedirs(output_dir)
    except OSError:
        if not isdir(output_dir):
            raise


def write_options_to_log(log, opts):
    """
    Writes user's input options to log file

    :param log: log filename
    :param opts: options
    """


    logfile = open(join(opts.output, log),"w+")
    logfile_header = "#Karenina Simulation Logfile\n"
    logfile.write(logfile_header)

    logfile.write("Output folder: %s\n" %(str(opts.output)))
    logfile.write("Treatment names: " + (str(opts.treatment_names)) + "\n")
    n_individuals_line = "Number of individuals: %s\n"\
    %(str(opts.n_individuals))
    logfile.write(n_individuals_line)
    logfile.write("Number of timepoints: " + (str(opts.n_timepoints)) + "\n")
    logfile.write("Perturbation timepoint: " +
    (str(opts.perturbation_timepoint)) + "\n")
    logfile.write("Perturbation duration: " +
    (str(opts.perturbation_duration)) + "\n")
    logfile.write("Interindividual variation: " +
    (str(opts.interindividual_variation)) + "\n")
    logfile.write("Delta: " + (str(opts.delta)) + "\n")
    logfile.write("Lambda: " + (str(opts.L)) + "\n")
    logfile.write("Fixed starting position: " + (str(opts.fixed_start_pos)) +
    "\n")

    logfile.close()



def parse_perturbation_file(pert_file_path, perturbation_timepoint,perturbation_duration):
    """
    Return a list of perturbations
    infile -- a .tsv file describing one perturbation per line
    assume input file is correctly formatted (no warnings if not)

    NOTE: each pertubation should be in the format:

    set_xyz_lambda_low = {"start":opts.perturbation_timepoint,
    "end":opts.perturbation_timepoint + opts.perturbation_duration,
    "params":{"lambda":0.005}, "update_mode":"replace", "axes":["x","y","z"]}

    :param pert_file_path: perturbation file path
    :param perturbation_timepoint: timepoint to apply perturbation
    :param perturbation_duration: duration of perturbation
    :return: perturbation list parsed from pert file contents
    """

    perturbations_list = []

    if (pert_file_path != None):
        df = pd.read_csv(pert_file_path, sep = "\t")

        headers_list = list(df)

        for index, row in df.iterrows():

            a_perturbation = {"start":perturbation_timepoint,\
            "end":perturbation_timepoint + perturbation_duration}

            required_headers_checker = {"params" : False, "values" : False,
            "update_mode" : False, "axes" : False}

            for header in headers_list:

                header_lowercase = header.lower()

                if header_lowercase in ("parameter", "parameters", "param",\
                "params"):
                    required_headers_checker["params"] = True
                    params = row[header].split(",")

                elif header_lowercase in ("value", "values", "val", "vals"):
                    required_headers_checker["values"] = True
                    values = str(row[header]).split(",")

                elif header_lowercase in ("update_mode", "update_modes",\
                "update mode", "update modes"):
                    required_headers_checker["update_mode"] = True
                    update_mode = row[header]

                elif header_lowercase in ("axes", "axis"):
                    required_headers_checker["axes"] = True
                    axes = row[header].split(",")

                else:
                    raise ValueError("Could not identify header name in " + \
                    "perturbations file")

            missing_headers_error_message = ""
            for each_checker in required_headers_checker:
                if required_headers_checker[each_checker] == False:
                    missing_headers_error_message += each_checker + " "
            if missing_headers_error_message != "":
                missing_headers_error_message = "Missing the following " +\
                "header(s): " + missing_headers_error_message
                raise ValueError(missing_headers_error_message)

            if len(params) != len(values):
                raise ValueError("Number of parameters does not match the " + \
                "number of values")

            a_perturbation["params"] = {}
            for idx, single_param in enumerate(params):
                a_perturbation["params"][single_param] = float(values[idx])
            a_perturbation["update_mode"] = update_mode
            a_perturbation["axes"] = axes

            perturbations_list.append(a_perturbation)

    else:
        set_xyz_lambda_zero = {"start":perturbation_timepoint,\
        "end":perturbation_timepoint + perturbation_duration,\
        "params":{"lambda":0.000},"update_mode":"replace","axes":["x","y","z"]}

        perturbations_list.append(set_xyz_lambda_zero)

    return perturbations_list


def main():

    parser = make_option_parser()

    opts, args = parser.parse_args()

    if opts.output is None:
        parser.print_help()
        exit()


    write_options_to_log("log.txt", opts)

    verbose = opts.verbose
    #Check timepoints
    check_perturbation_timepoint(opts.perturbation_timepoint,opts.n_timepoints)
    #Set the base parameters for microbiome change over time
    #in unperturbed individuals.
    individual_base_params = {"lambda":opts.L,"delta":opts.delta,\
      "interindividual_variation":opts.interindividual_variation}
    if opts.fixed_start_pos:
        try:
            x,y,z = map(float,opts.fixed_start_pos.split(","))
            individual_base_params['x']=x
            individual_base_params['y']=y
            individual_base_params['z']=z

        except:
            print ("Supplied value for fixed start position after parsing:",opts.fixed_start_pos)
            raise ValueError('Problem with --fixed_start_pos. Got %s Please supply x,y,z values in the range (-1,1) separated by commas and enclosed in quotes. Example: "0.1,-0.2,0.3"'% opts.fixed_start_pos)

    #Set up the treatments to be applied

    perturbations = parse_perturbation_file(opts.pert_file_path,\
    opts.perturbation_timepoint, opts.perturbation_duration)

    treatments = [[], perturbations]
    treatment_names = opts.treatment_names.split(",")
    if verbose:
        print("Raw number of individuals from user:",opts.n_individuals)
        print("n_individuals",opts.n_individuals.split(','))
    n_individuals = list(map(int,opts.n_individuals.split(",")))
    if verbose:
        print ("**Experiment Design**")
        print ("treatments:",treatment_names)
        print ("n_individuals:",n_individuals)
        print ("interindividual_variation",opts.interindividual_variation)
        print ("treatment_effects:",treatments)
        print ("individual_base_params:",individual_base_params)
    experiment = Experiment(treatment_names,n_individuals,opts.n_timepoints,\
        individual_base_params,treatments,opts.interindividual_variation, verbose)
    experiment.simulate_timesteps(0,opts.n_timepoints, verbose)
    experiment.write_to_movie_file(opts.output, verbose)

if __name__ == "__main__":
    main()
