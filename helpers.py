import json
from datetime import datetime, timedelta

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


from plot_tno import plot_df

def load_config(path):
    """
    Loads the configuration from json file and returns it as dict
    :param path: path to json file
    :return: dict representation of json
    """
    with open(path, 'r') as f:
        config = json.load(f)
    return config


def save_results(param, results, outpath):
    with h5py.File(outpath, 'w') as hf:
        parameters = hf.create_group('parameters')
        res = hf.create_group('results')
        for p in param:
            if p == 'time':
                res[p] = param[p]
            else:
                # check for np
                try:
                    if np.min(param[p]) == np.max(param[p]):
                        parameters[p] = np.min(param[p])
                    else:
                        out = np.array(param[p])
                        if len(out.shape) == 2:
                            parameters[p] = param[p][0]
                        else:
                            parameters[p] = param[p]
                except:
                    pass


        names = ['x', 'inf',  'ft', 'Rt', 'Rtd', 'Rtratio', 'pu', 'pv', 'pb', 'infpeak', 'hosp', 'hosppeak']
        results = [np.array([r[a] for r in results]) for a in range(len(results[0]))]
        #for i, name in enumerate(names[:-1]):
        for i, name in enumerate(names):
            res[name] = results[i]


def load_results(h5path):
    parameters = {}
    results = {}
    with h5py.File(h5path, 'r') as hf:
        param_keys = hf['parameters'].keys()
        for key in param_keys:
            parameters[key] = hf['parameters'][key][()]

        result_keys = hf['results'].keys()
        for key in result_keys:
            results[key] = hf['results'][key][()]

    return parameters, results


def visualize_result_summary(h5path, figpath, plotname, plotcolor, plottitle, x=None, plothistograms=True,
                             plotscale =1.0,  tabulate=False, plotdata="", plotshift=0):
    parameters, results = load_results(h5path)

    try:
        assert len(x) > 0
    except:
        x = results['time'][0]



    nr_realizations = results['infpeak'].size
    varying_param = [key for key in parameters if np.size(parameters[key]) == nr_realizations]
    temp = results[plotname]
    plotscale2 =1.0
    if (plotscale>=0):
        plotscale2 =plotscale
    mean_t = np.mean(temp, axis=0)*plotscale2
    median_t = np.percentile(temp, 50, axis=0)*plotscale2
    p5 = np.percentile(temp, 5, axis=0)*plotscale2
    p95 = np.percentile(temp, 95, axis=0)*plotscale2

    fig = plt.figure(figsize=(20, 5), constrained_layout=False)

    outer_grid = fig.add_gridspec(1, 4)

    row_col_count = int(np.ceil(np.sqrt(len(varying_param))))
    if (plothistograms):
        inner_grid = outer_grid[0].subgridspec(row_col_count, row_col_count, wspace=0.5, hspace=0.5)
        for i, par in enumerate(varying_param):
            ax = fig.add_subplot(inner_grid[i])
            ax.hist(parameters[par])
            ax.set_xlabel(par)
            fig.add_subplot(ax)
    """   
    for i, par in enumerate(varying_param):
        plt.subplot(row_col_count, row_col_count, i + 1)
        plt.scatter(parameters[par], hospeak)
        plt.xlabel(par)
        plt.ylabel('hospeak')
    """
    ax = fig.add_subplot(outer_grid[1:])
    x = x + timedelta(days=plotshift)
    ax.plot(x, mean_t, c='k', label='Mean')
    ax.plot(x, median_t, c='k', ls=':', label='P50')
    ax.fill_between(x, p5, p95, facecolor=plotcolor, label='90% conf')
    ax.set_xlabel('Date/days')
    ax.set_ylabel(plottitle)
    if (plotdata!=""):
        data = pd.read_csv(plotdata, delimiter='\t')
        data.day = pd.to_datetime(data.day, format="%d-%m-%Y")
        x_obs = data.day.values
        y_obs = data.cases.values*1.0
        if (plotscale>=0):
            y_obs/= y_obs[plotshift]
        plt.scatter(x_obs, y_obs , c='k', label='data', marker='o', s=8)


    #if plotname == 'Rt':
    #    ax.set_ylim(1, 1.5)
    plt.title(plottitle)
    #ax.set_ylim(0,10)
    plt.grid(True)
    ax.legend()
    fig.add_subplot(ax)

    if (tabulate):
        print(plotname)
        print ( 'day, value')
        for i, val in enumerate(x):
            print ( i, mean_t[i])

    sumfigpath = figpath + plotname + '.png'
    plt.savefig(sumfigpath, dpi=300)
    plt.close()


def visualize_cross_plots(h5path, figpath, skey, scale=1.0):
    parameters, results = load_results(h5path)

    nr_realizations =  results[skey].size
    hospeak = results[skey]
    varying_param = [key for key in parameters if np.size(parameters[key]) == nr_realizations]
    row_col_count = int(np.ceil(np.sqrt(len(varying_param))))

    for i, par in enumerate(varying_param):

        plt.subplot(row_col_count, row_col_count, i + 1)
        plt.scatter(parameters[par], hospeak*scale, alpha=0.5)
        plt.xlabel(par)
        if (scale==1.0):
            plt.ylabel('daily hospitalized (index)')
        else:
            plt.ylabel('daily hospitalized')

    plt.tight_layout()
    plt.savefig(figpath, dpi=300)
    plt.close()

def read_df(h5path, skey):
    parameters, results = load_results(h5path)
    nr_realizations =  results[skey].size
    varying_param = [key for key in parameters if np.size(parameters[key]) == nr_realizations]
    hospeak = results[skey]
    data = { skey : hospeak}
    df = pd.DataFrame(data)
    for i, par in enumerate(varying_param):
        data = parameters[par]
        df[par] = data
    return df

def visualize_crossplot_df(h5path, skey):
    df = read_df(h5path, skey)
    parameters, results = load_results(h5path)
    nr_realizations = results[skey].size
    varying_param = [key for key in parameters if np.size(parameters[key]) == nr_realizations]
    for i, par in enumerate(varying_param):
        plot_df(df, par, skey, percentile_window = 10,percentile_values=(90,70), percentile_fill=(True, True), percentile_colors=("mistyrose","powderblue"),
                percentile_alpha=0.8)

def visualize_cross_plots_plevels(h5path, figpath, skey, plevels, plevelcolor, scale=1.0):
    parameters, results = load_results(h5path)

    nr_realizations =  results[skey].size
    hospeak = results[skey]
    varying_param = [key for key in parameters if np.size(parameters[key]) == nr_realizations]
    row_col_count = int(np.ceil(np.sqrt(len(varying_param))))

    ymax = np.max(hospeak)
    a = 1 - 0.2*np.log(nr_realizations)/np.log(10)
    nbin = (int) (nr_realizations/1000.0)
    for i, par in enumerate(varying_param):
        #ax = plt.add_subplot(row_col_count, row_col_count, i + 1)
        ax = plt.subplot(row_col_count, row_col_count, i + 1)
        plt.ylim(0,ymax*scale)
        data = parameters[par]
        s = np.min(data)
        e = np.max(data)
        bins = np.linspace(s, e, nbin+1)
        binindex = np.digitize(data, bins)
        binhospeaks = [np.sort(hospeak[binindex==i]) for i in range(1, len(bins))]

        for iplevel, plevel in enumerate(plevels):
            bins2 = bins[:-1]
            xbin = np.empty([bins2.size * 2])
            p5 = np.empty([bins2.size *2])
            p95 =  np.empty([bins2.size *2])

            for i, pp in enumerate(binhospeaks):
                xbin[i*2] = bins[i]
                xbin[i*2 +1] = bins[i+1]
                p5[i*2] = np.percentile(binhospeaks[i], 0)
                p5[i*2+1] = p5[i*2]
                p95[i*2] = np.percentile(binhospeaks[i], 100- plevel)
                p95[i * 2 + 1] = p95[i * 2]
            ax.fill_between(xbin, p5*scale, p95*scale, facecolor=plevelcolor[iplevel], alpha=0.8, label='90% conf')
        plt.scatter(parameters[par], hospeak*scale, alpha=a, s=0.5)
        plt.plot()
        plt.grid(True)
        plt.xlabel(par)
        if (scale==1.0):
            plt.ylabel('daily (index=1)')
        else:
            plt.ylabel('daily hospitalized')


    #plt.legend(loc = 'lower right')
    plt.tight_layout()
    if (scale == 1.0):
        sumfigpath = figpath + skey + 'index.png'
    else:
        sumfigpath = figpath + skey + '.png'
    plt.savefig(sumfigpath, dpi=300)
    plt.close()
