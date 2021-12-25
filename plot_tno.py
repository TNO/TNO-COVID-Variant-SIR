import os.path as path

import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np


SAVE_PLOTS = 0
DISPLAY_PLOTS = 1
SAVE_AND_DISPLAY_PLOTS = 2


def get_p_names(perc):
    sname1 = "P" + str(perc)
    sname2 = "P" + str(100 - perc)
    return sname1, sname2






def get_percentiles(df, xname, yname, percentile_window=3, percentile_values=(50, 90)):
        """

        :param df: dataframe, it  assumes the data has been sorted on xname by df.sort_values(xname)
        :param xname: x coordinate field
        :param yname: y coordinate field
        :param percentile_window:   size of rolling window to pot percentile values, -1 means no rolling window results to plot
        :param percentile_values : percentile values (0..100.0) to plot with a moving window in xname
               it constructs two lines from
                 df[xname].rolling(window = percentile_window, center = True).quantile(percentile_value/100)
                 df[xname].rolling(window = percentile_window, center = True).quantile(1-percentile_value/100)
        :return copy of df with
        """
        if percentile_window < 0:
            return
        else:
            if percentile_window == 0:
                percentile_window = np.ceil(np.size(df[xname]) / 25)
                percentile_window = max(3, percentile_window)
            for i, perc in enumerate(percentile_values):
                sname1, sname2 = get_p_names(perc)
                df[sname1] = df[yname].rolling(window=percentile_window, center=True, min_periods=percentile_window).quantile(
                    perc * 1.0 / 100.0)
                df[sname2] = df[yname].rolling(window=percentile_window, center=True, min_periods=percentile_window).quantile(
                    (100 - perc) * 1.0 / 100.0)
            #  df[sname1] = df[sname1].rolling(window=percentile_window, center=True).mean()
            #  df[sname2] = df[sname2].rolling(window=percentile_window, center=True).mean()




def plot_df(df_input: pd.DataFrame, xname: str, yname: str, catcol=None, size=None, alpha=1.0, input_figsize=(10, 10),
                ylog=False, newplot=True, show=True, input_plotsizes=(40, 400), percentile_window=-1, percentile_values=None,
                percentile_fill=None, percentile_colors=None, percentile_alpha=0.8):
        """ plots a pandas dataframe in an x,y plot. The symbols are circles and their colors and sizes can be
        controlled by color and size fields. For future versions plotting of this type can be done with seaborn

        :param df_input: dataframe
        :param xname: x coordinate field
        :param yname: y coordinate field
        :param catcol: field to denote colors, this can be a non-number field, such that it generates categories
        :param size: field to denote size
        :param alpha:  opacity (0 is fully transparant, 1 is fully opaque)
        :param input_figsize:  size in cm of plot x-axis  and y-axis (float)
        :param ylog:  show logaritmic scale (boolean)
        :param newplot:
        :param show
        :param input_plotsizes: minimum and maximum symbol area size
        :param percentile_window:   size of rolling window to plot percentile values, -1 means no rolling window results to plot
        :param percentile_values : percentile values (0..100.0) to plot with a moving window in xname
               it constructs two lines from
                 df[xname].rolling(window = percentile_window, center = True).quantile(percentile_value/100)
                 df[xname].rolling(window = percentile_window, center = True).quantile(1-percentile_value/100)
        :param  percentile_fill boolean to plot lines or fill
        :param  percentile_colors to be used for lines and fill
        :param  percentile_alpha fixed alpha value for fill

        """

        df = df_input.copy()
        df = df.sort_values(xname)
        get_percentiles(df, xname, yname, percentile_window, percentile_values)


        # x = np.asarray(df[xname])
        # y = np.asarray(df[yname])

        x = df[xname]
        y = df[yname]
        if np.size(x) == 0:
            return

        if newplot:
            fig, axes = plt.subplots(1, 1, figsize=input_figsize)
        else:
            plt.gcf()
            axes = plt.gca()

        plt.xlabel(xname)
        plt.ylabel(yname)
        if ylog:
            plt.yscale("log")

        # start with percentile plots if any
        if percentile_window > 0:
            for i, perc in enumerate(percentile_values):
                if percentile_fill[i]:
                    sname1, sname2 = get_p_names(perc)
                    slabel = sname1
                    axes.fill_between(x, df[sname1], df[sname2], facecolor=percentile_colors[i], alpha=percentile_alpha,
                                      label=slabel)

        plotsize_max = input_plotsizes[1]
        plotsizes_scale = plotsize_max / np.log10(np.size(x))
        plotsizes_scalemin = input_plotsizes[0] / input_plotsizes[1]
        plotsizes_scalemax = (input_plotsizes[1] - input_plotsizes[0]) / input_plotsizes[1]
        if size is not None:
            sizes = np.asarray(df[size])
            sizes2 = (sizes - min(sizes)) / (max(sizes) - min(sizes))
            df['plotsizes'] = plotsizes_scale * (plotsizes_scalemin + plotsizes_scalemax * sizes2)
        else:
            df['plotsizes'] = input_plotsizes[0]

        ncolors = 6
        categories = False
        groups = None
        colors = None
        if catcol is not None:
            try:
                dummy = df[catcol] * 1.0
                colors = np.asarray(df[catcol])
            except:
                categories = True
                groups = df.groupby(catcol)
        else:
            colors = None

        scatter = None
        if categories:
            for name, group in groups:
                xg = np.asarray(group[xname])
                yg = np.asarray(group[yname])
                plotsizesg = np.asarray(group['plotsizes'])
                if np.size(x) != 0:
                    scatter = plt.scatter(xg, yg, s=plotsizesg, alpha=alpha, label=name)
        else:
            plotsizes = np.asarray(df['plotsizes'])
            scatter = plt.scatter(x, y, c=colors, s=plotsizes, alpha=alpha, cmap="jet")

        if percentile_window > 0:
            # plotlines
            for i, perc in enumerate(percentile_values):
                if not (percentile_fill[i]):
                    sname1, sname2 = get_p_names(perc)
                    plt.plot(x, df[sname1], color=percentile_colors[i], linewidth=1)
                    plt.plot(x, df[sname2], color=percentile_colors[i], linewidth=1)


        if categories:
            legend1 = plt.legend(loc="lower right", title=catcol)
            axes.add_artist(legend1)
        else:
            if catcol is not None:
                legend1 = axes.legend(*scatter.legend_elements(num=ncolors),
                                      loc="lower right", title=catcol)
                axes.add_artist(legend1)

        if size is not None:
            kw = dict(prop="sizes", num=6, color=scatter.cmap(0.7), fmt=" {x:.2f}",
                      func=lambda s: (((s / plotsizes_scale) - plotsizes_scalemin) / plotsizes_scalemax) * (
                              max(sizes) - min(sizes)) + min(sizes))
            axes.legend(*scatter.legend_elements(**kw), loc="upper right", title=size)

        plt.grid()
        if show:
            plt.show()
        else:
            return plt
