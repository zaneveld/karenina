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

from experiment import Experiment
import visualization
from optparse import OptionParser
from optparse import OptionGroup
from os.path import join,isdir,realpath,dirname
from os import makedirs


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
      default = join(dirname(dirname(realpath(__file__))),'data',\
      'perturbations','xyz_lambda_zero.csv'),\
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

    parser.add_option_group(optional_options)

    return parser


def check_perturbation_timepoint(perturbation_timepoint,n_timepoints):
    """Raise ValueError if perturbation_timepoint is < 0 or >n_timepoints"""
    if perturbation_timepoint and perturbation_timepoint >= n_timepoints:
        raise ValueError("Perturbation timepoint must be before the last timepoint")
    if perturbation_timepoint < 0:
        raise ValueError("Perturbation timepoint must be positive")


def ensure_exists(output_dir):
    """Ensure that output_dir exists"""
    try:
        makedirs(output_dir)
    except OSError:
        if not isdir(output_dir):
            raise


def write_options_to_log(log, opts):
    """Writes user's input options to log file"""

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


def parse_perturbation_file(opts):
    """Return a list of perturbations
    infile -- a .csv file describing one perturbation per line
    assume input file is correctly formatted (no warnings if not)

    NOTE: each pertubation should be in the format:
    set_xyz_lambda_low =
       {"start":opts.perturbation_timepoint,\
       "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"lambda":0.005},\
      "update_mode":"replace",\
      "axes":["x","y","z"]}
    """
    perturb_list = []
    if (opts.pert_file_path != None):
        input_file = open(opts.pert_file_path, "r")

        for line in input_file:
            print (line)
            temp_list = []

            for word in line.split('\t'):
                temp_list.append(word)

            temp_dict = {"start":opts.perturbation_timepoint,\
              "end":opts.perturbation_timepoint + opts.perturbation_duration}

            for temp_list_index in range(len(temp_list)):
                if (temp_list[temp_list_index] == 'params'):
                    index = temp_list_index + 1
                    temp_dict['params'] = {}
                    while (temp_list[index] != 'update_mode'):
                        temp_dict['params'][temp_list[index]] =\
                        float(temp_list[index + 1])
                        index += 2
                elif (temp_list[temp_list_index] == 'update_mode'):
                    temp_dict["update_mode"] = temp_list[temp_list_index + 1]
                elif (temp_list[temp_list_index] == 'axes'):
                    index = temp_list_index + 1;
                    temp_dict["axes"] = []
                    while (index != len(temp_list) and temp_list[index] != ''
                    and temp_list[index] != '\n'):
                        split = temp_list[index].split('\n')
                        temp_dict["axes"].append(split[0])
                        index += 1
            perturb_list.append(temp_dict)
        input_file.close()

    else:
        set_xyz_lambda_zero = {"start":opts.perturbation_timepoint,\
        "end":opts.perturbation_timepoint + opts.perturbation_duration,\
        "params":{"lambda":0.000},"update_mode":"replace","axes":["x","y","z"]}

        perturb_list.append(set_xyz_lambda_zero)

    return perturb_list


def main():

    parser = make_option_parser()
    opts, args = parser.parse_args()
    print (opts)

    write_options_to_log("log.txt", opts)

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

    perturbations = parse_perturbation_file(opts)

    treatments = [[], perturbations]
    treatment_names = opts.treatment_names.split(",")
    print("Raw number of individuals from user:",opts.n_individuals)
    print("n_individuals.split(',')",opts.n_individuals.split(','))
    n_individuals = list(map(int,opts.n_individuals.split(",")))

    print ("**Experiment Design**")
    print ("treatments:",treatment_names)
    print ("n_individuals:",n_individuals)
    print ("interindividual_variation",opts.interindividual_variation)
    print ("treatment_effects:",treatments)
    print ("individual_base_params:",individual_base_params)
    experiment = Experiment(treatment_names,n_individuals,opts.n_timepoints,\
        individual_base_params,treatments,opts.interindividual_variation)
    experiment.simulate_timesteps(0,opts.n_timepoints)
    experiment.writeToMovieFile(opts.output)

if __name__ == "__main__":
    main()
