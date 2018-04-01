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

from spatial_ornstein_uhlenbeck_individual import Individual
from spatial_ornstein_uhlenbeck_perturbation import Perturbation
from os.path import join,isdir,realpath,dirname
from numpy import array
from copy import copy


def get_timeseries_data(individuals,axes=["x","y","z"]):
    results = []
    for i,curr_subject in enumerate(individuals):
        result = []
        for j,axis in enumerate(axes):
            entry = curr_subject.MovementProcesses[axis].History
            result.append(entry)
        results.append(array(result))
    return results


def save_simulation_figure(individuals, output_folder,n_individuals,n_timepoints,perturbation_timepoint):
    """Save a .pdf image of the simulated PCoA plot"""

    individual_colors = {"healthy":"orange","perturbed":"magenta"}
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(1,1,1) # one row, one column, first plot
    #ax.set_axis_bgcolor('black') deprecated
    ax.set_facecolor('black')
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

    data = get_timeseries_data(individuals)
    colors = [i.BaseParams["color"] for i in individuals]
    print("Individual colors:",colors)
    print ("Movie raw data:",data)
    # NOTE: Can't pass empty arrays into 3d version of plot()
    linestyle = '-'
    pointstyle = 'o' #cheat to use lines to represent points
    lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1],linestyle,\
      c=colors[i],alpha=0.20)[0] for i,dat in enumerate(data)]

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
        #ax.set_axis_bgcolor('black') deprecated
        ax.set_facecolor('black')
        fig.patch.set_facecolor('black')
        dull_red = (0.50,0,0,1)
        ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 1))
        ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 1))
        ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 1))
        #Set 3d background grid color to a dull red
        new_grid_params = ax.w_xaxis._axinfo['grid']
        new_grid_params.update({'color': dull_red, 'linewidth':1.0})
        print(new_grid_params)
        ax.w_xaxis._axinfo.update({'grid' : new_grid_params})
        ax.w_yaxis._axinfo.update({'grid' : new_grid_params})
        ax.w_zaxis._axinfo.update({'grid' : new_grid_params})

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

        self.TreatmentNames = [t for t in treatment_names]
        self.Treatments = [{"treatment_name":name} for name in self.TreatmentNames]
        self.BaseParams = individual_base_params
        self.NIndividuals = [n for n in n_individuals]
        self.NTimepoints = n_timepoints
        #Check that a few parameters are valid
        print ("treatment_names:",self.TreatmentNames)
        print ("n_individuals:",self.NIndividuals)
        print ("treatment_params:",treatment_params)
        self.check_n_timepoints_is_int(n_timepoints)
        self.check_variable_specified_per_treatment(self.NIndividuals)
        self.check_variable_specified_per_treatment(treatment_params)


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

            print(treatment)
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
        print([x for x in v])
        print(self.TreatmentNames)
        if len([x for x in v]) != len(self.TreatmentNames):
            raise ValueError('Must specify a list of n_individuals equal in length to the number of treatments. Note that n_individuals must be enclosed in quotes e.g. -n "35,35"  ')

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

