import argparse
import os
import sys
import time
import numpy as np
from dateutil import rrule, parser
import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # To ensure proper importing below
from helpers import load_config, save_results, visualize_result_summary, visualize_cross_plots, \
    visualize_cross_plots_plevels, visualize_crossplot_df
from plot_tno import plot_df
from voccore import run_models




def main(args):
    """
    Run VOC analysis

    :param args: All arguments passed to the main script. First argument has to be path to the configuration .json
                 Other arguments are not currently in use
    :return: None
    """
    start_time = time.time()

    # Parse the configuration file and combine with standard configuration
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("configfile_path", type=str, help='The full path to the .json config file')
    options = None
    if args is not None:
        options, arguments = arg_parser.parse_known_args(args=args)
    config = load_config(options.configfile_path)

    param, results = run_models(config)

    print('Saving results and visualizing...')
    h5path = os.path.join(config['base_dir'],config['run_name'])+'.h5'
    sumfigpath = os.path.join(config['base_dir'],config['run_name'])
    crsfigpath = os.path.join(config['base_dir'], config['run_name']) + '_crossplot.png'
    save_results(param, results, h5path)


    datesx = list(rrule.rrule(rrule.DAILY, dtstart=parser.parse(config['startdate']), until=parser.parse(config['enddate'])))
    x = np.array(datesx)

    plotnames = config['plot']['plotnames']
    plotcolors= config['plot']['plotcolors']
    plottitles = config['plot']['plottitles']
    plothistograms = config['plot']['plothistograms']
    plotscales = config['plot']['scales']
    plottimeshift= config['plot']['timeshift']
    plotdata = config['plot']['data']

    for i, plotname in enumerate(plotnames):
        visualize_result_summary(h5path, sumfigpath, plotname, plotcolors[i], plottitles[i], x, plothistograms,plotscales[i],
                                 tabulate= (plotname=='Rt'),plotdata=plotdata[i], plotshift=plottimeshift[i])

    visualize_cross_plots(h5path, crsfigpath, "infpeak")
    plevels = [10, 25, 50]
    plevelcolors = ["red", "orange","yellow"]

    visualize_cross_plots_plevels(h5path,sumfigpath, "infpeak",plevels, plevelcolors, scale=1.0 )
    visualize_cross_plots_plevels(h5path, sumfigpath, "hosppeak", plevels, plevelcolors, scale=1.0)
    #visualize_crossplot_df(h5path, "hosppeak")


    print('Total runtime: {:.2f} minutes'.format((time.time() - start_time) / 60.))


if __name__ == '__main__':
    main(sys.argv[1:])
