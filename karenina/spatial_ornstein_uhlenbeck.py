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



from scipy.stats import norm
from sys import argv
from random import random,randint
#from cogent.util.option_parsing import parse_command_line_parameters, make_option
from optparse import OptionParser
from optparse import OptionGroup
from ast import literal_eval
import numpy as np
import matplotlib.pyplot as plt
from os.path import join,isdir
from os import makedirs
from numpy import array
from copy import copy

"""
script_info = {}
script_info['brief_description'] = "This script simulates microbial community change over time in PCoA-space using Brownian Motion or Ornstein-Uhlenbeck models."
script_info['script_description'] = "This script simulates microbiome change over time using Ornstein-Uhlenbeck (OU) models.  These are similar to Brownian motion models, with the exception that they include reversion to a mean"
script_info['script_usage'] = [
                               ("","Simulate microbiomes using default parameters .", "%prog -o ./simulation_results")
                                ]
script_info['output_description']= "Output is a tab-delimited data table and figures."

script_info['required_options'] = [
 make_option('-o','--output',type="new_filepath",help='the output folder for the simulation results')
]

script_info['optional_options'] = [\
    make_option('--treatment_names',default="control,destabilizing_treatment",type="string",help="Comma seperated list of treatment named [default:%default]"),
    make_option('-n','--n_individuals',default="35,35",type="string",help='Comma-separated number of individuals to simulate per treatment. [default: %default]'),
    make_option('-t','--n_timepoints',default=10,type="int",help='Number of timepoints to simulate. (one number, which is the same for all treatments) [default: %default]'),
    make_option('-p','--perturbation_timepoint',default=5,type="int",help='Timepoint at which to apply a perturbation. Must be less than --n_timepoints [default: %default]'),
    make_option('-d','--perturbation_duration',default=100,type="int",help='Duration that the perturbation lasts. [default: %default]'),
    make_option('--interindividual_variation',default=0.01,type="float",help='Starting variability between individuals. [default: %default]'),
    make_option('--delta',default=0.25,type="float",help='Starting delta parameter for Brownian motion and Ornstein-Uhlenbeck processes. A higher number indicates more variability over time. [default: %default]'),
    make_option('-l','--L',default=0.20,type="float",help='Starting lambda parameter for Ornstein-Uhlenbeck processes. A higher number indicates a greater tendancy to revert to the mean value. [default: %default]'),
    make_option('--fixed_start_pos',default=None,type="string",help='Starting x,y,z position for all points, as comma separated floating point values, e.g. 0.0,0.1,0.2. If not supplied, starting positions will be randomized based on the interindividual_variation parameter [default: %default]')
    ]

script_info['version'] = __version__
"""

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

    """ OptionGroup gives error: invalid option type??
    required_options.add_option('-o','--output',type="new_filepath",
    help='the output folder for the simulation results')
    """

    parser.add_option_group(required_options)


    optional_options = OptionGroup(parser, "Optional options")

    optional_options.add_option("-i", "--input_file",\
    action="store", type="string",\
    help = "Input file for analysis")

    optional_options.add_option('--treatment_names',\
    default="control,destabilizing_treatment",type="string",\
    help="Comma seperated list of treatment named [default:%default]")

    optional_options.add_option('-n','--n_individuals',\
    default="35,35",type="string",\
    help='Comma-separated number of individuals to simulate per treatment. ' +
    '[default: %default]')

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

class Process(object):

    #Constructor to Process class
    #class variables don't have to be declared first?
    """Represents a 1d process in a Euclidean space"""
    def __init__(self,start_coord, motion = "Ornstein-Uhlenbeck",\
        history = None,params={"L":0.20,"delta":0.25}):
        """
        start_coords - float starting coordinate for the particle
        """
        if history is None:
            history = []
        self.StartCoord = start_coord
        self.Coord = start_coord
        self.History = history
        self.History.append(start_coord)
        self.Params = params #hash table?
        self.ProcessType = motion #string
        self.Perturbations = [] #this is an empty list

    def update(self,dt):
        curr_params = copy(self.Params) #deep or shallow copy?
        for p in self.Perturbations:
            curr_params = p.updateParams(curr_params) #function defined below
        if self.ProcessType == "Brownian":
            self.bm_update(dt,delta=curr_params["delta"])

            #if it's not Brownian then it must be OU? then else should be sufficient?
        elif self.ProcessType == "Ornstein-Uhlenbeck": #does the \ mean including the next line?
            self.ou_update(dt,mu=curr_params["mu"],\
            delta = curr_params["delta"],\
            L=curr_params["lambda"])

    #what is bm?
    def bm_change(self,dt,delta):
        change =  norm.rvs(loc=0,size=1,scale=delta**2*dt) #what does this do?
        return change

    def bm_update(self,dt,delta):
        curr_coord = self.Coord
        self.History.append(curr_coord)
        change = self.bm_change(dt,delta)
        self.Coord = curr_coord + change


    def ou_change(self,dt,mu,\
        L,delta):
        """
        The Ornstein Uhlenbeck process is modelled as:

        ds = lambda * (mu - s) * dt + dW

        ds -- change in our process from the last timepoint
        L -- lambda, the speed of reversion to mean
        (NOTE: lambda is a reserved keyword in Python so I use L)
        mu -- mean position
        s -- current position
        dt -- how much time has elapsed since last update
        dW -- the Weiner Process (basic Brownian motion)

        This says we update as usual for Brownian motion,
        but add in a term that reverts us to some
        mean position (mu) over time (dt) at some speed (lambda)

        """

        dW = self.bm_change(dt=dt,delta=delta)
        ds = L * (mu - self.Coord) * dt + dW
        return ds

    def ou_update(self,dt,mu,\
        L,delta,min_bound=-1.0,max_bound=1.0):
        curr_coord = self.Coord
        self.History.append(self.Coord)
        change = self.ou_change(dt=dt,mu=mu,L=L,delta=delta)
        self.Coord = curr_coord + change
        if min_bound is not None:
            self.Coord = max(self.Coord,min_bound)
        if max_bound is not None:
            self.Coord = min(self.Coord,max_bound)


class Perturbation(object):
    def __init__(self,start,end,params,update_mode="replace",axes=["x","y","z"]):
        """Alter a simulation to impose a press disturbance,
           shifting the microbiome torwards a new configuration

           start --inclusive timepoint to start perturbation.  Note that this is read at the
             Experiment level, not by underlying Process objects.

           end -- inclusive timepoint to end perturbation. Note that this is read at the
             Experiment level, not by underlying Process objects.

           params -- dict of parameter values altered by the disturbance

           mode -- how the perturbation updates parameter values.
                'replace' -- replace old value with new one
                'add' -- add new value to old one
                'multiply' -- multiply the two values

           axes -- axes to which the perturbation applies.  Like Start and End this is a
            'dumb' value, read externally by the Experiment object
        """
        self.Start = start
        self.End = end
        self.Params = params
        self.UpdateMode = update_mode
        self.Axes = axes

    def isActive(self,t):
        return self.Start <= t <= self.End

    def updateParams(self,params):
        if self.UpdateMode == "replace":
            update_f = self.updateByReplacement
        elif self.UpdateMode == "add":
            update_f = self.updateByAddition
        elif self.UpdateMode == "multiply":
            update_f = self.updateByMultiplication
        else:
            raise ValueError("Invalid update mode for perturbation: %s" %self.UpdateMode)

        #If we set new_params = params
        #we modify the input dict, which
        #isn't what we want
        new_params = copy(params)
        for k,v in self.Params.iteritems():
            new_params[k] = update_f(params[k],v)

        return new_params

    def updateByReplacement(self,curr_param, perturbation_param):
        return perturbation_param

    def updateByAddition(self,curr_param, perturbation_param):
        return curr_param + perturbation_param

    def updateByMultiplication(self,curr_param, perturbation_param):
        return curr_param * perturbation_param

class Individual(object):
    def __init__(self,subject_id,coords=["x","y","z"],metadata={},params={},interindividual_variation=0.01):
        self.SubjectId = subject_id
        self.Metadata = metadata
        self.MovementProcesses = {}
        self.BaseParams = params
        for c in coords:
            #print "STARTING PROCESS for Axis:",c

            #simulate minor interindividual variation
            if c in params.keys():
                start_coord = params[c]
                start_mu = params[c]
            else:
                start_coord =  (random()-0.50) *\
                     2.0 * ( interindividual_variation)
                start_mu = start_coord
            #print "START COORD %s: %f" %(c,start_coord)

            #For OU processes, assume the process reverts
            # to its starting coordinate
            curr_params = copy(self.BaseParams)
            curr_params["mu"] = start_coord
            print ("start_coord:",start_coord)
            print ("curr_params['mu']",curr_params['mu'])
            self.MovementProcesses[c] = Process(start_coord = start_coord,params=curr_params,\
              motion = "Ornstein-Uhlenbeck")

    def applyPerturbation(self,perturbation):
        """Apply a perturbation to the appropriate axes"""
        for axis in perturbation.Axes:
            self.applyPerturbationToAxis(axis,perturbation)

    def removePerturbation(self,perturbation):
        for axis in perturbation.Axes:
            self.removePerturbationFromAxis(axis,perturbation)

    def removePerturbationFromAxis(self,axis,perturbation):
        """Remove a perturbation from one or more Process objects"""
        self.MovementProcesses[axis].Perturbations.remove(perturbation)

    def applyPerturbationToAxis(self,axis,perturbation):
        """Apply a perturbation to a Processes objects"""

        self.MovementProcesses[axis].Perturbations.append(perturbation)

    def check_identity(self):
        for coord_name,movement_process in self.MovementProcesses.iteritems():
            for coord_name2,movement_process2 in self.MovementProcesses.iteritems():
                if coord_name == coord_name2:
                    continue
                if movement_process.History == movement_process2.History:
                    #print "History1:",movement_process.History
                    #print "History2:",movement_process2.History
                    return True
        return False

    def simulate_movement(self,n_timepoints,params=None):
        if not params:
            params = self.BaseParams

        for c in self.MovementProcesses.keys():
            for t in range(n_timepoints):
                mu = self.MovementProcesses[c].StartCoord
                self.MovementProcesses[c].update(dt=1.0)

    def get_data(self,n_timepoints):
        result = []
        coords = self.MovementProcesses.keys()
        for t in range(n_timepoints):
            timepoint = [self.MovementProcesses[c].History[t] for c in coords]
            curr_sampleid = "S%s_t%i"%(self.SubjectId,t)
            curr_data = [curr_sampleid]
            curr_data.extend(map(float,timepoint))
            result.append(curr_data)
        return result

def save_simulation_figure(individuals, output_folder,n_individuals,n_timepoints,perturbation_timepoint):
    """Save a .pdf image of the simulated PCoA plot"""

    individual_colors = {"healthy":"orange","perturbed":"magenta"}


    # Create a Figure object.
    fig = plt.figure(figsize=(5, 4))
    # Create an Axes object.
    ax = fig.add_subplot(1,1,1) # one row, one column, first plot
    ax.set_axis_bgcolor('black')
    fig.patch.set_facecolor('black')
    ax.set_title('Simulated Microbiome Destabilization\n (n=%i; %i individuals, %i timepoints)' %(n_timepoints * n_individuals,\
      n_individuals,n_timepoints), color='w', fontsize=12)
    ax.set_xlabel(r'PC1 (simulated)')
    ax.set_ylabel('PC2 (simulated)')
    ax.set_xlim((-1.0,1.0))
    ax.set_ylim((-1.0,1.0))
    plt.figtext(0.55,0.75,r"BM: $\frac{\delta x}{\delta t} = \sigma W_t$",color="orange")
    plt.figtext(0.55,0.65,r"OU: $\frac{\delta x}{\delta t} = \sigma W_t + (\Theta - x)\lambda  $",color="magenta")
    ax.spines['bottom'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    series = []
    for i,curr_subject in enumerate(individuals):
        #Plot pre-perturbation timepoints
        xs = curr_subject.MovementProcesses["x"].History[:perturbation_timepoint]
        ys = curr_subject.MovementProcesses["y"].History[:perturbation_timepoint]
        curr_color =individual_colors["healthy"]
        series_handle = ax.scatter(xs,ys,c=curr_color,s=36,edgecolor=None,alpha=0.5)
        series.append(series_handle)

        #Plot post-perturbation timepoints
        xs = curr_subject.MovementProcesses["x"].History[perturbation_timepoint:]
        ys = curr_subject.MovementProcesses["y"].History[perturbation_timepoint:]
        curr_color =individual_colors["perturbed"]
        series_handle = ax.scatter(xs,ys,c=curr_color,s=36,edgecolor=None,alpha=0.5)
        series.append(series_handle)


    fig_filename = join(output_folder,"simulation_%i_hosts_%i_timepoints.pdf" %(n_individuals,n_timepoints))
    fig.savefig(fig_filename, facecolor=fig.get_facecolor(), edgecolor='none',bbox_inches='tight')


def get_timeseries_data(individuals,n_timepoints,start=0,end=-1,axes=["x","y","z"]):
    results = []
    for i,curr_subject in enumerate(individuals):
        result = []
        for j,axis in enumerate(axes):
            entry = curr_subject.MovementProcesses[axis].History
            result.append(entry)
        results.append(array(result))
    return results

def update_3d_plot(end_t,timeseries_data,ax,lines,points=None,start_t=0):

    for line,data in zip(lines,timeseries_data):
        line.set_data(data[0:2,start_t:end_t])
        #z pos can't be set with set_data
        line.set_3d_properties(data[2,start_t:end_t])

    if points:
         for point,data in zip(points,timeseries_data):
            point.set_data(data[0:2,end_t-1:end_t])
            #z pos can't be set with set_data
            point.set_3d_properties(data[2,end_t-1:end_t])
    rotation_speed = 0.5
    ax.view_init(30, rotation_speed * end_t)

def save_simulation_movie(individuals, output_folder,\
     n_individuals,n_timepoints,\
    black_background=True):
    """Save an .ffmpg move of the simulated community change"""

    #TODO: standardize these and put them up above

    import numpy as np
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d.axes3d as p3
    import matplotlib.animation as animation


    #The code for writing animation files is essentially identical to the
    #matplotlib tutorial here: http://matplotlib.org/examples/animation/basic_example_writer.html
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist=str(__author__)), bitrate=1800)

    # Attaching 3D axis to the figure
    fig = plt.figure()
    ax = p3.Axes3D(fig)

    data = get_timeseries_data(individuals,0, n_timepoints)
    colors = [i.BaseParams["color"] for i in individuals]
    print("Individual colors:",colors)
    print ("Movie raw data:",data)
    # NOTE: Can't pass empty arrays into 3d version of plot()
    linestyle = '-'
    pointstyle = 'o' #cheat to use lines to represent points
    lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1],linestyle,\
      c=colors[i],alpha=0.20)[0] for i,dat in enumerate(data)]

    pointstyle = 'o'
    points = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1],pointstyle,\
      c=colors[i],alpha=1.0)[0] for i,dat in enumerate(data)]

    # Setting the axes properties
    ax.set_xlim3d([-1.0, 1.0])
    ax.set_xlabel('PC1')

    ax.set_ylim3d([-1.0, 1.0])
    ax.set_ylabel('PC2')

    ax.set_zlim3d([-1.0, 1.0])
    ax.set_zlabel('PC3')

    ax.set_title('Simulation Results')
    if black_background:
        ax.set_axis_bgcolor('black')
        fig.patch.set_facecolor('black')
        dull_red = (0.50,0,0,1)
        ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 1))
        ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 1))
        ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 1))
        #Set 3d background grid color to a dull red
        ax.w_xaxis._axinfo.update({'grid' : {'color': dull_red}})
        ax.w_yaxis._axinfo.update({'grid' : {'color': dull_red}})
        ax.w_zaxis._axinfo.update({'grid' : {'color': dull_red}})

        ax.spines['bottom'].set_color(dull_red)
        ax.spines['left'].set_color(dull_red)
        ax.xaxis.label.set_color(dull_red)
        ax.yaxis.label.set_color(dull_red)
        ax.zaxis.label.set_color(dull_red)
        ax.tick_params(axis='x', colors=dull_red)
        ax.tick_params(axis='y', colors=dull_red)
        ax.tick_params(axis='z', colors=dull_red)

    # Creating the Animation object
    line_ani = animation.FuncAnimation(fig, update_3d_plot, n_timepoints, fargs=(data,ax,lines,points),\
      interval=100, blit=False)
    line_ani.save(join(output_folder,'simulation_video.mp4'), writer=writer)
    #plt.show()




class Experiment(object):
    """This class has responsibility for simulating an experimental design for simulation.

    A fixed number of individuals are simulated across experimental conditions called 'treatments'
    Each treatment can have different numbers of individuals.
    All treatment must have the same number of timepoints.

    Each treatment can be associated with the imposition of one or more Perturbations.
    Each perturbation is inserted or removed from all individuals at a fixed time-point.


    """
    def __init__(self,treatment_names,n_individuals,n_timepoints,\
        individual_base_params,treatment_params,interindividual_variation):
        """Set up an experiment with multiple treatments

        Parameters
        ----------
        treatments -- a list of treatment names. e.g. ['Control','Temperature Stress']
        n_individuals -- a list of the number of individuals in each treatment. e.g. [10,4]
        n_timepoints -- an integer number of timepoints
        individual_base_params -- a default dict with the base parameters (pre-perturbation) to be
          applied to each individual.
        treatment_params -- a list of lists of perturbations to apply to each treatment. The form is a bit
          complex to allow for more powerful specifications of experimental design.  See note below.

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

        interindividual_varation -- the amount of starting variation between individuals
        """

        self.TreatmentNames = treatment_names
        self.Treatments = [{"treatment_name":name} for name in self.TreatmentNames]
        self.BaseParams = individual_base_params
        self.NIndividuals = n_individuals
        self.NTimepoints = n_timepoints
        #Check that a few parameters are valid
        print ("treatment_names:",treatment_names)
        print ("n_individuals:",n_individuals)
        print ("treatment_params:",treatment_params)
        self.check_n_timepoints_is_int(n_timepoints)
        self.check_variable_specified_per_treatment(n_individuals)
        self.check_variable_specified_per_treatment(treatment_params)


        for i,n in enumerate(n_individuals):
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

            for i in range(treatment["n_individuals"]):

                curr_subject_id = "%s_%i" %(treatment["treatment_name"],i)
                curr_subject = Individual(subject_id = curr_subject_id,
                  params = params,\
                  metadata={"treatment":treatment["treatment_name"]},\
                  interindividual_variation=interindividual_variation)
                individuals.append(curr_subject)
            treatment["individuals"] = individuals


        #Set up the treatment parameters
        for treatment_idx,treatment in enumerate(self.Treatments):
            treatment["perturbations"] = []
            treatment["active_perturbations"] = []
            raw_perturbation_info = treatment_params[treatment_idx]
            #We should have a dict with start, end, and parms for each perturbation
            print ("raw_perturbation_info:",raw_perturbation_info)
            for p in raw_perturbation_info:
                print ("params:",p)
                curr_perturbation = Perturbation(p["start"],p["end"],p["params"],p["update_mode"],p["axes"])
                treatment["perturbations"].append(curr_perturbation)

        #Set up a place to hold data on the experiment outcome

        coords = ["x","y","z"]
        headers = "\t".join(["SampleID"]+coords)+"\n"
        self.Data = [headers]

    def run(self):
        "Run the experiment, simulating timesteps"
        save_simulation_figure(individuals,opts.output,n_individuals,n_timepoints,perturbation_timepoint)
        save_simulation_movie(individuals, opts.output,n_individuals,n_timepoints,perturbation_timepoint)

    def check_variable_specified_per_treatment(self,v):
        """Raise a ValueError if v is not the same length as the number of treatments"""
        if len(v) != len(self.TreatmentNames):
            raise ValueError("Must specify a list of n_individuals equal in length to the number of treatments")

    def check_n_timepoints_is_int(self,n_timepoints):
        """Raise a ValueError if n_timepoints can't be cast as an int"""
        try:
            n_timepoints = int(n_timepoints)
        except:
            raise ValueError("n_timepoints must be a single integer that applies to all experiments (not a list per treatment for example).")

    def simulate_timesteps(self,t_start,t_end):
        """Simulate multiple timesteps"""
        for t in range(t_start,t_end):
            print ("Simulating timestep: %i" %t)
            self.simulate_timestep(t)

    def simulate_timestep(self,t):
        """Simulate timestep t of the experiemnt

        Approach:
        First apply any perturbations that should be active but aren't yet
        Second, simulate the timestep
        Third, remove any perturbations that should be off

        NOTE: this means that a perturbation that starts and ends at t=1 will be active at t=1
        That is, the start and end times are inclusive.
        """

        for treatment in self.Treatments:
            for perturbation in treatment["perturbations"]:
                apply_to_individuals = False
                remove_from_individuals = False
                #Record whether to activate perturbation
                if perturbation.isActive(t) and\
                    perturbation not in treatment["active_perturbations"]:
                    apply_to_individuals = True
                    treatment["active_perturbations"].append(perturbation)

                #Record whether to deactivate perturbation
                if perturbation in treatment["active_perturbations"] and\
                    not perturbation.isActive(t):
                    remove_from_individuals = True
                    treatment["active_perturbations"].remove(perturbation)

                #Apply new perturbations and remove old ones
                for curr_subject in treatment["individuals"]:
                    if apply_to_individuals:
                        curr_subject.applyPerturbation(perturbation)
                    if remove_from_individuals:
                        curr_subject.removePerturbation(perturbation)

            #With perturbations in place and updated,
            #simulate a timestep for each individual in treatment
            for curr_subject in treatment["individuals"]:
                    #Simulate the timestep
                    curr_subject.simulate_movement(1)
                    curr_data = curr_subject.get_data(1)
                    self.Data.append("\t".join(map(str,curr_data))+"\n")

    def writeToMovieFile(self,output_folder):
        """Write an MPG movie to output folder"""
        individuals = []
        for treatment in self.Treatments:
            for curr_subject in treatment["individuals"]:
                individuals.append(curr_subject)
        print ("individuals:",individuals)
        save_simulation_movie(individuals, output_folder,\
             len(individuals),self.NTimepoints,\
             black_background=True)


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

def main():

    """
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    ensure_exists(opts.output)

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
            print (fixed_start_pos)
            raise ValueError("Problem with --fixed_start_pos. Got %s Please supply tx,y,z values in the range (-1,1) separated by commas. Example: 0.1,-0.2,0.3"% fixed_start_pos)
    """
    #Set up the treatments to be applied

    #TODO: parameterize this by parsing a parameter file
    set_xyz_lambda_low = {"start":opts.perturbation_timepoint,\
       "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"lambda":0.005},"update_mode":"replace","axes":["x","y","z"]}

    set_yz_lambda_medium = {"start":opts.perturbation_timepoint,\
       "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"lambda":0.08},"update_mode":"replace","axes":["z","y"]}

    set_y_lambda_medium = {"start":opts.perturbation_timepoint,\
       "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"lambda":0.16},"update_mode":"replace","axes":["y"]}

    set_xyz_lambda_zero = {"start":opts.perturbation_timepoint,\
       "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"lambda":0.000},"update_mode":"replace","axes":["x","y","z"]}

    set_x_mu_low = {"start":opts.perturbation_timepoint,\
       "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"mu":-0.8},"update_mode":"replace","axes":["x"]}

    double_xyz_delta  = {"start":opts.perturbation_timepoint,\
      "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"delta":2.0},"update_mode":"multiply","axes":["x","y","z"]}

    double_z_delta  = {"start":opts.perturbation_timepoint,\
      "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"delta":2.0},"update_mode":"multiply","axes":["z"]}

    set_xyz_mu_low = {"start":opts.perturbation_timepoint,\
      "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"mu":-0.8},"update_mode":"replace","axes":["x","y","z"]}

    set_xyz_mu_high = {"start":opts.perturbation_timepoint,\
      "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"mu":0.8},"update_mode":"replace","axes":["x","y","z"]}

    add_x_mu_high = {"start":opts.perturbation_timepoint,\
      "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"mu":0.8},"update_mode":"add","axes":["x"]}

    set_x_mu_high = {"start":opts.perturbation_timepoint,\
      "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"mu":0.8},"update_mode":"replace","axes":["x"]}

    set_x_lambda_small = {"start":opts.perturbation_timepoint,\
      "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"lambda":0.008,"mu":-0.08},"update_mode":"replace","axes":["x"]}

    set_x_lambda_medium = {"start":opts.perturbation_timepoint,\
      "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"lambda":0.016,"mu":-0.08},"update_mode":"replace","axes":["x"]}

    set_x_lambda_high = {"start":opts.perturbation_timepoint,\
      "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"lambda":0.064},"update_mode":"replace","axes":["x"]}

    set_yz_lambda_high = {"start":opts.perturbation_timepoint,\
      "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"lambda":0.064},"update_mode":"replace","axes":["z","y"]}

    set_y_lambda_high = {"start":opts.perturbation_timepoint,\
      "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"lambda":0.064},"update_mode":"replace","axes":["y"]}

    set_z_lambda_zero = {"start":opts.perturbation_timepoint,\
      "end":opts.perturbation_timepoint + opts.perturbation_duration,\
      "params":{"lambda":0.0},"update_mode":"replace","axes":["z"]}

    #TODO let user choose
    #treatments = [[],[set_x_lambda_high]]
    #treatments = [[],[double_xyz_delta]]
    #treatments = [[],[set_xyz_lambda_low]]
    #treatments = [[],[add_x_mu_high,double_z_delta]]
    #treatments = [[],[set_xyz_mu_high]]
    #treatments = [[],[double_xyz_delta]]
    treatments = [[],[set_xyz_lambda_zero]]
    treatment_names = opts.treatment_names.split(",")
    n_individuals = map(int,opts.n_individuals.split(","))

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
