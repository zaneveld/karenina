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
from karenina.fit_timeseries import parse_pcoa,parse_metadata
from optparse import OptionParser
from optparse import OptionGroup
from os.path import join
from numpy import array
from scipy.spatial import distance
import os
import re
import random
import csv


def make_option_parser():
    """Return an optparse OptionParser object"""

    parser = OptionParser(usage = "%prog -o ./simulation_results",
        description="This script fits microbiome change over time to Ornstein-Uhlenbeck (OU) models."+
                    "Demo fitting an OU model using default parameters.",
        version= __version__)

    required_options = OptionGroup(parser, "Required options")

    required_options.add_option('-o','--output', type="string",
                                help='the output folder for the simulation results')

    parser.add_option_group(required_options)

    optional_options = OptionGroup(parser, "Optional options")

    optional_options.add_option('--pcoa_qza', type="string", default=None,
                                help='The input PCoA qza file')
    optional_options.add_option('--metadata', type="string", default=None,
                                help='The input metadata tsv file, if not defined, '
                                     'metadata will be extracted from PCoA qza')

    # Individual identifier will either accept a single column, or optionally combine two columns of categorical features.
    optional_options.add_option('--individual', type="string", default=None,
                                help='The host subject ID column identifier, separate by comma to combine TWO columns')
    optional_options.add_option('--timepoint', type="string", default=None,
                                help='The timepoint column identifier')
    optional_options.add_option('--treatment', type="string", default=None,
                                help='The treatment column identifier')

    optional_options.add_option('-v', '--verbose', action="store_true", dest="verbose", default=False,
                                help='-v, allows for verbose output' +
                                     ' [default: %default]')

    parser.add_option_group(optional_options)

    return parser

def get_timeseries_data(individuals,axes=["x","y","z"]):
    """
    Provides timeseries data from list of individuals.

    :param individuals: array of individual objects
    :param axes: [x,y,z] axes to return
    :return: array of timeseries data
    """
    results = []
    for i,curr_subject in enumerate(individuals):
        result = []
        for j,axis in enumerate(axes):
            entry = curr_subject.MovementProcesses[axis].History
            result.append(entry)
        results.append(array(result))
    return results


def save_simulation_figure(individuals, output_folder,n_individuals,n_timepoints,perturbation_timepoint):
    """
    Save a .pdf image of the simulated PCoA plot

    :param individuals: array of individuals
    :param output_folder: output filepath
    :param n_individuals: number of individuals
    :param n_timepoints: number of timepoints
    :param perturbation_timepoint: timepoint of perturbation application
    """

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


def save_simulation_data(data, ids, output):
    """
    Saves simulation output data in PCoA format

    :param data: data to save
    :param ids: Sample_IDs
    :param output: output filepath
    """
    with open(output+"ordination.txt","w") as outfile:
        # Need to calculate eigenvalues
        # outfile.write("Eigvals\t" + str(len(data)) + "\n\n")
        outfile.write("Eigvals\t0" + "\n\n")

        # Need to calculate propEx
        # outfile.write("Proportion explained\t" + str(len(data)) + "\n\n")
        outfile.write("Proportion explained\t0"+ "\n\n")

        outfile.write("Species\t0\t0\n\n")
        outfile.write("Site\t"+str(len(data)*len(data[0][0]))+"\t3\n")

        # Need to separate pc1,2,3 and assign unique identifiers based on hash and timepoint.
        dm = {}
        for j,row in enumerate(data):
            identifier = ids[j]
            for i in range(len(row[0])):
                outfile.write(str(identifier)+"_t"+str(i)+"\t"+str(row[0][i])+"\t"+str(row[1][i])+"\t"+str(row[2][i])+"\n")
                dm.update({str(identifier)+"."+str(i):[row[0][i],row[1][i],row[2][i]]})

        outfile.write("\n")
        outfile.write("Biplot\t0\t0\n\n")
        outfile.write("Site constraints\t0\t0\n")
    outfile.close()

    # Distance matrix (euclidean)
    dm_0 = []
    dm_0.append("")
    distance_matrix = []
    for key in dm.keys():
        dm_0.append(key)
    distance_matrix.append(dm_0)
    for key in dm.keys():
        dm_1 = []
        dm_1.append(key)
        for key1 in dm.keys():
            dm_1.append(str(distance.euclidean(dm[key],dm[key1])))
        distance_matrix.append(dm_1)

    with open(output+"euclidean.txt","w") as outfile:
        for row in distance_matrix:
            for item in row:
                outfile.write(str(item)+"\t")
            outfile.write("\n")
    outfile.close()

    #Mapping file
    md_0 = ["#SampleID","Subject","Treatment","Timepoint"]
    md_1 = ["#q2:types","categorical","categorical","numeric"]
    md = []
    for id in ids:
        for i in range(len(data[0][0])):
            md.append([id+"_t"+str(i),id,''.join([k for k in id if not k.isdigit()])[:-1],i])
    metadata = [md_0,md_1]
    for row in md:
        metadata.append(row)
    with open(output+"metadata.tsv","w") as outfile:
        for row in metadata:
            for i,item in enumerate(row):
                if i < len(row)-1:
                    outfile.write(str(item)+"\t")
                if i==len(row)-1:
                    outfile.write(str(item))
            outfile.write("\n")
    outfile.close()


def save_simulation_movie(individuals, output_folder,\
     n_individuals,n_timepoints,black_background=True, data_o = False, verbose=False):
    """
    Save an .ffmpg move of the simulated community change

    :param individuals: array of individuals to visualize
    :param output_folder: output directory filepath
    :param n_individuals: number of individuals
    :param n_timepoints: number of timepoints
    :param black_background: T/F, default = True
    :param verbose: verbose output, default = False
    """

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
    if verbose:
        print("Individual colors:",colors)
        print("Movie raw data:",data)

    ids = []
    for item in individuals:
        ids.append(item.SubjectId)
    if data_o:
        save_simulation_data(data, ids, output_folder)
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
        if verbose:
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
    """
    Updates visualization 3d plot

    :param end_t: end timepoint
    :param timeseries_data: data from timeseries
    :param ax: visualization ax
    :param lines: lines of data
    :param points: values to update
    :param start_t: start timepoint (0)
    """
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


def main():
    parser = make_option_parser()
    opts, args = parser.parse_args()
    verbose = opts.verbose
    if verbose:
        print(opts)

    if opts.output is None:
        raise IOError('Output must be defined.')

    if not os.path.exists(opts.output):
        os.makedirs(opts.output)

    tx = opts.treatment
    output = opts.output
    n_timepoints = opts.timepoint

    ind = []

    site, metadata = parse_pcoa(opts.pcoa_qza, opts.individual, opts.timepoint, opts.treatment, opts.metadata)
    df = parse_metadata(metadata, opts.individual, opts.timepoint, opts.treatment, site)

    colors = ['fuchsia', 'cyan', 'darkorange', 'blue', 'yellow']
    tx = opts.treatment
    treatments = df[tx].unique()
    while len(colors) < len(treatments):
        colors.append('lightgray')

    for i,row in enumerate(df.iterrows()):
        curr_subject_id = "%s_%i" % (df[opts.individual], i)
        j=0
        while row[1][3] != treatments[j]:
            j+=1
        color = colors[j]
        params = {'lambda': 0.2, 'delta': 0.25, 'interindividual_variation': 0.01}
        params['color'] = color
        row_tx = row[1][3]
        curr_subject = Individual(subject_id=curr_subject_id,
                                  params=params, \
                                  metadata={opts.treatment: row_tx}, \
                                  interindividual_variation=.01, verbose=verbose)
        ind.append(curr_subject)


    #save_simulation_figure(individuals=ind, output_folder=output, n_timepoints=50, perturbation_timepoint=25, n_individuals=50)
    save_simulation_movie(individuals=ind,
                          output_folder=output,
                          n_timepoints=len(df[opts.timepoint].unique()),
                          n_individuals=len(df[opts.individual].unique()),
                          verbose=verbose)

if __name__ == "__main__":
    main()