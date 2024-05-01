import sys
import os
import re
import numpy as np
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import Divider, Size
from statsmodels.nonparametric.smoothers_lowess import lowess
import argparse


bar_width = 0.8
capthick = 1
capsize = 2
x_sep = 1
font_size = 16
tick_size = 14
legend_size = 10
rotation = 30
orangered = colors.to_rgba('orangered')
dodgerblue = colors.to_rgba('dodgerblue')
gray = colors.to_rgba('gray')
gray = tuple([x if i < 3 else 0.4 for i, x in enumerate(gray)])
plt.rcParams.update({'font.size': font_size, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size})
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"


for parcel_set in ('evlab', 'PDD', 'RH', 'PDDanat'):
    ylims = {'': {1: (-0.05, 1.1), 2: (-0.05, 0.35)}}
    if parcel_set == 'evlab':
        fROIs = [
            'LIFGorb',
            'LIFG',
            'LMFG',
            'LAntTemp',
            'LPostTemp',
            'LAngG'
        ]
    elif parcel_set.startswith('PDD'):
        if parcel_set == 'PDDanat':
            ylims['_tight'] = {1: (-0.05, 0.5), 2: (-0.02, 0.19)}
        fROIs = [
            'LIFGorb',
            'LIFGtri',
            'LTP',
            'LaSTS',
            'LpSTS',
            'LTPJ',
        ]
    else:
        fROIs = [
            'RIFGorb',
            'RIFG',
            'RMFG',
            'RAntTemp',
            'RPostTemp',
            'RAngG'
        ]

    try:
        df1 = pd.read_csv('contrasts/%s_nlength1_contrasts.csv' % parcel_set)
        df2 = pd.read_csv('contrasts/%s_nlength2_contrasts.csv' % parcel_set)
    except FileNotFoundError:
        sys.stderr.write('Contrast files not found. Run `python -m nlength.contrasts` before running this script.\n')
        sys.stderr.flush()
        exit()

    # Main result

    cmap = plt.get_cmap('terrain')
    hatches = [
        None,
        '\\\\\\\\',
        '||||',
        '----',
        '....',
        'xxxx',
        '++++'
    ]
    colors = [cmap(i/len(hatches)) for i in range(len(hatches))]
    bar_width = 1./ 6 * 0.5


    fig = plt.figure(figsize=((4, 4.5)))
    h = [Size.Fixed(0.7), Size.Fixed(3.3)]
    v = [Size.Fixed(0.5), Size.Fixed(3.7)]
    divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
    ax = fig.add_axes(
        divider.get_position(),
        axes_locator=divider.new_locator(nx=1, ny=1)
    )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.tick_params(labelleft='on', labelbottom='off')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    contrast = 'C_v_J'


    for i, fROI in enumerate(['all'] + fROIs):
        r = i * bar_width

        if fROI == 'all':
            df = df2[contrast]
        else:
            df = df2[df2.fROI == fROI][contrast]
        if contrast == 'C_v_J':
            df = df
        mean = float(df.mean())
        err = float(df.sem())

        ax.bar(
            r,
            mean,
            # color='w',
            color=colors[i],
            edgecolor='none',
            # edgecolor=(0.6, 0.6, 0.6),
            width=bar_width,
            label='Overall' if fROI == 'all' else fROI,
            # hatch=hatches[i],
            linewidth=2
        )

        ax.errorbar(
            r,
            mean,
            yerr=err,
            fmt='none',
            ecolor=colors[i],
            # ecolor=(0.8, 0.8, 0.8),
            capsize=0,
            capthick=2,
            linewidth=2
        )

    ax.set_xticks([])
    ax.set_xlim(-0.25, 0.75)
    ax.axhline(y=0, color='k', lw=0.5)

    for ylim_key in ylims:
        ylim = ylims[ylim_key][1]
        ax.set_ylim(ylim)
        plt.savefig('plots/%s_overall1%s.png' % (parcel_set, ylim_key))
    plt.close('all')


    fig = plt.figure(figsize=((14, 4.5)))
    h = [Size.Fixed(0.7), Size.Fixed(13.3)]
    v = [Size.Fixed(0.5), Size.Fixed(3.7)]
    divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
    ax = fig.add_axes(
        divider.get_position(),
        axes_locator=divider.new_locator(nx=1, ny=1)
    )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.tick_params(labelleft='on', labelbottom='off')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    contrasts = ['NLenC', 'NLenJ', 'NLenC_v_NLenJ']
    contrasts_renamed = [''] * (len(contrasts) + 1)

    r_base = np.arange(len(contrasts_renamed))

    for i, fROI in enumerate(['all'] + fROIs):
        r = r_base + i * bar_width
        means = []
        errs = []

        if fROI == 'all':
            df = df1[['NLenC']]
        else:
            df = df1[df1.fROI == fROI][['NLenC']]
        means.append(float(df.mean()))
        errs.append(float(df.sem()))

        for contrast in contrasts:
            if fROI == 'all':
                df = df2[contrast]
            else:
                df = df2[df2.fROI == fROI][contrast]
            if contrast == 'C_v_J':
                df = df
            means.append(float(df.mean()))
            errs.append(float(df.sem()))

        ax.bar(
            r,
            means,
            # color='w',
            color=colors[i],
            edgecolor='none',
            # edgecolor=(0.6, 0.6, 0.6),
            width=bar_width,
            label='Overall' if fROI == 'all' else fROI,
            # hatch=hatches[i],
            linewidth=2
        )

        ax.errorbar(
            r,
            means,
            yerr=errs,
            fmt='none',
            ecolor=colors[i],
            # ecolor=(0.8, 0.8, 0.8),
            capsize=0,
            capthick=2,
            linewidth=2
        )

    ax.legend(bbox_to_anchor=(1,1.1), ncol=4, fontsize=20)

    ax.set_xticks([])
    ax.axhline(y=0, color='k', lw=0.5)

    for ylim_key in ylims:
        ylim = ylims[ylim_key][2]
        ax.set_ylim(ylim)
        plt.savefig('plots/%s_overall2%s.png' % (parcel_set, ylim_key))
    plt.close('all')