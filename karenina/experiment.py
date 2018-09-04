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

from karenina.individual import Individual
from karenina.perturbation import Perturbation
import karenina.visualization as visualization
from copy import copy



class Experiment(object):
    """
    This class has responsibility for simulating an experimental design for simulation.

    A fixed number of individuals are simulated across experimental conditions called 'treatments'
    Each treatment can have different numbers of individuals.
    All treatment must have the same number of timepoints.

    Each treatment can be associated with the imposition of one or more Perturbations.
    Each perturbation is inserted or removed from all individuals at a fixed time-point.


    """
    def __init__(self,treatment_names,n_individuals,n_timepoints,\
        individual_base_params,treatment_params,interindividual_variation, verbose):
        """
        Set up an experiment with multiple treatments

        Specifying perturbations:
        Perturbations are specified as a dict with parameters matching a 'Perturbation' object: start,end,params,update_mode representing the
        starting and ending timesteps of the perturbation, the parameters of that perturbation, and how the specified parameters affect the
        individuals pre-existing parameters (e.g. by replacement or addition).  The
        parameters of the perturbation are specified in a nested dict by axis. Mode can be set to "add", "multiply" or "replace"

        So for example:

        set_lambda_low_treatment = {"start":10, "end:"25,"params":{"x":{"lambda":0.005},"mode":"replace"}}
        treatment_names = ["control","destablizing_treatment"]
        treatments = [[],[set_lambda_low_treatment]]

        -- this will specify a two treatment experiment, in which the 'destabilizing treatment

        :param treatment_names: a list of treatment names. e.g. ['Control','Temperature Stress']
        :param n_individuals: a list of the number of individuals in each treatment. e.g. [10,4]
        :param n_timepoints: an integer number of timepoints
        :param individual_base_params: a default dict with the base parameters (pre-perturbation) to be
          applied to each individual.
        :param treatment_params: a list of lists of perturbations to apply to each treatment. The form is a bit
          complex to allow for more powerful specifications of experimental design.
        :param interindividual_variation: the amount of starting variation between individuals
        :param verbose: verbose output
        """
        self.verbose = verbose
        self.TreatmentNames = [t for t in treatment_names]
        self.Treatments = [{"treatment_name":name} for name in self.TreatmentNames]
        self.BaseParams = individual_base_params
        self.NIndividuals = [n for n in n_individuals]
        self.NTimepoints = n_timepoints
        #Check that a few parameters are valid
        if self.verbose:
            print("treatment_names:",self.TreatmentNames)
            print("n_individuals:",self.NIndividuals)
            print("treatment_params:",treatment_params)
        self.check_n_timepoints_is_int(n_timepoints)
        self.check_variable_specified_per_treatment(self.NIndividuals, self.verbose)
        self.check_variable_specified_per_treatment(treatment_params, self.verbose)


        #for i,n in enumerate(n_individuals):
        for i,n in enumerate(self.NIndividuals):
            self.Treatments[i]["n_individuals"] = n

        self.NTimepoints = n_timepoints

        colors = ['fuchsia','cyan','darkorange','blue','yellow']
        #Set up the experimental subjects
        for treatment_idx,treatment in enumerate(self.Treatments):

            #Set a color for each treatment
            individuals = []
            params = copy(self.BaseParams)

            if treatment_idx < len(colors):
                params['color'] = colors[treatment_idx]
            else:
                params['color'] = 'lightgray'

            if self.verbose:
                print(treatment)
            for i in range(treatment["n_individuals"]):
                curr_subject_id = "%s_%i" %(treatment["treatment_name"],i)
                curr_subject = Individual(subject_id = curr_subject_id,
                  params = params,\
                  metadata={"treatment":treatment["treatment_name"]},
                  interindividual_variation=interindividual_variation, verbose=self.verbose)
                individuals.append(curr_subject)
            treatment["individuals"] = individuals


        #Set up the treatment parameters
        for treatment_idx,treatment in enumerate(self.Treatments):
            treatment["perturbations"] = []
            treatment["active_perturbations"] = []
            raw_perturbation_info = treatment_params[treatment_idx]
            #We should have a dict with start, end, and parms for each perturbation
            if self.verbose:
                print("raw_perturbation_info:",raw_perturbation_info)
            for p in raw_perturbation_info:
                #
                if self.verbose:
                    print("params:",p)
                curr_perturbation = Perturbation(p["start"],p["end"],p["params"],p["update_mode"],p["axes"])
                treatment["perturbations"].append(curr_perturbation)

        #Set up a place to hold data on the experiment outcome

        coords = ["x","y","z"]
        headers = "\t".join(["SampleID"]+coords)+"\n"
        self.Data = [headers]

    def run(self):
        """
        Run the experiment, simulating timesteps; saves the simulation figure and movie
        """
        visualization.save_simulation_figure(individuals,opts.output,n_individuals,n_timepoints,perturbation_timepoint)
        visualization.save_simulation_movie(individuals, opts.output,n_individuals,n_timepoints,perturbation_timepoint, verbose)

    def check_variable_specified_per_treatment(self,v,verbose=False):
        """
        Raise a ValueError if v is not the same length as the number of treatments

        :param v: variable
        :param verbose: verbose output, default = False
        """
        if verbose:
            print([x for x in v])
            print(self.TreatmentNames)
        if len([x for x in v]) != len(self.TreatmentNames):
            raise ValueError('Must specify a list of n_individuals equal in length to the number of treatments. Note that n_individuals must be enclosed in quotes e.g. -n "35,35"  ')

    def check_n_timepoints_is_int(self,n_timepoints):
        """
        Raise a ValueError if n_timepoints can't be cast as an int

        :param n_timepoints: number of timepoints
        """
        try:
            self.n_timepoints = int(n_timepoints)
        except:
            raise ValueError("n_timepoints must be a single integer that applies to all experiments (not a list per treatment for example).")

    def simulate_timesteps(self,t_start,t_end,verbose=False):
        """
        Simulate multiple timesteps

        :param t_start: start timepoint
        :param t_end: end timepoint
        :param verbose: verbose output, default = False
        """
        for t in range(t_start,t_end):
            if verbose:
                print ("Simulating timestep: %i" %t)
            self.simulate_timestep(t)



    def simulate_timestep(self,t):
        """
        Simulate timestep t of the experiemnt

        Approach:

        1. First apply any perturbations that should be active but aren't yet
        2. Second, simulate the timestep
        3. Third, remove any perturbations that should be off

        NOTE: this means that a perturbation that starts and ends at t=1 will be active at t=1
        That is, the start and end times are inclusive.

        :param t: timestep
        """

        for treatment in self.Treatments:
            for perturbation in treatment["perturbations"]:
                apply_to_individuals = False
                remove_from_individuals = False
                #Record whether to activate perturbation
                if perturbation.is_active(t) and\
                    perturbation not in treatment["active_perturbations"]:
                    apply_to_individuals = True
                    treatment["active_perturbations"].append(perturbation)

                #Record whether to deactivate perturbation
                if perturbation in treatment["active_perturbations"] and\
                    not perturbation.is_active(t):
                    remove_from_individuals = True
                    treatment["active_perturbations"].remove(perturbation)

                #Apply new perturbations and remove old ones
                for curr_subject in treatment["individuals"]:
                    if apply_to_individuals:
                        curr_subject.apply_perturbation(perturbation)
                    if remove_from_individuals:
                        curr_subject.remove_perturbation(perturbation)

            #With perturbations in place and updated,
            #simulate a timestep for each individual in treatment
            for curr_subject in treatment["individuals"]:
                #Simulate the timestep
                curr_subject.simulate_movement(1)
                curr_data = curr_subject.get_data(1)
                self.Data.append("\t".join(map(str,curr_data))+"\n")

    def write_to_movie_file(self,output_folder, verbose=False):
        """
        Write an MPG movie to output folder

        :param output_folder: output directory
        :param verbose: verbose output, default = False
        """
        individuals = []
        for treatment in self.Treatments:
            for curr_subject in treatment["individuals"]:
                individuals.append(curr_subject)
            #print("individuals:",individuals)

        for individual in individuals:
            individual.MovementProcesses["x"].History.pop(0)
            individual.MovementProcesses["y"].History.pop(0)
            individual.MovementProcesses["z"].History.pop(0)

        for item in individuals:
            print(vars(item))
        visualization.save_simulation_movie(individuals, output_folder,
                                            len(individuals),self.NTimepoints,
                                            black_background=True,
                                            data_o = True, verbose=self.verbose)
    def q2_data(self):
        """
        generate output data from object's self for Qiime2

        :return: data, ids
        """
        individuals = []
        for treatment in self.Treatments:
            for curr_subject in treatment["individuals"]:
                individuals.append(curr_subject)
            #print("individuals:",individuals)

        for individual in individuals:
            individual.MovementProcesses["x"].History.pop(0)
            individual.MovementProcesses["y"].History.pop(0)
            individual.MovementProcesses["z"].History.pop(0)

        data = visualization.get_timeseries_data(individuals)
        ids = []

        for item in individuals:
            ids.append(item.SubjectId)

        return data,ids