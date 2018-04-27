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

from os.path import join
from numpy import array

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
