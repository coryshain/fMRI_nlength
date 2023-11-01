import pickle
import os
import re
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import Divider, Size
import argparse

from nlength.contrasts import parcel_set_map, fROIs, networks
from nlength.tree import Tree


REGNUM = re.compile('(.+)_region([0-9]+)')
RUNNUM = re.compile('run([0-9]+)')
CLEN = re.compile('[^0-9]*([0-9]+)')


def z(x):
    return (x - x.mean()) / x.std()


def var(m, weights):
    w = np.array(weights)
    while len(w.shape) < 2:
        w = w[..., None]
    cov = m.cov_params()

    return np.dot(w.T, np.dot(cov, w))


def parse_path(path, kind='behavioral'):
    if kind.lower() == 'behavioral':
        p_ix = 3
        r_ix = 5
    elif kind.lower() == 'fmri':
        p_ix = 1
        r_ix = 2
    else:
        raise ValueError('Unrecognized path kind "%s".' % kind)

    path = path[:-4].strip().split('/')[-1]
    fields = path.split('_')
    subject = fields[p_ix]
    run = int(RUNNUM.match(fields[r_ix]).group(1))

    return subject, run


def cond2type(c):
    if c.startswith('cond'):
        out = 'C'
    elif '_jab' in c:
        out = 'JAB'
    elif c.endswith('nc'):
        out = 'NC'
    elif c.endswith('c'):
        out = 'C'
    else:
        out = c

    return out

def cond2len(c):
    match = CLEN.match(c)
    if match:
        return int(match.group(1))
    return 0


def cond2code(t, l):
    if t == 'C':
        if l == 12:
            out = 'A'
        elif l == 6:
            out = 'B'
        elif l == 4:
            out = 'C'
        elif l == 3:
            out = 'E'
        elif l == 2:
            out = 'G'
        elif l == 1:
            out = 'H'
        else:
            raise ValueError('Unrecognized length %d for condition %s.' % (l, t))
    elif t == 'NC':
        if l == 4:
            out = 'D'
        elif l == 3:
            out = 'F'
        else:
            raise ValueError('Unrecognized length %d for condition %s.' % (l, t))
    elif t == 'JAB':
        if l == 12:
            out = 'I'
        elif l == 4:
            out = 'J'
        elif l == 1:
            out = 'K'
        else:
            raise ValueError('Unrecognized length %d for condition %s.' % (l, t))
    else:
        out = 'FIX'

    return out


def get_docid(x):
    if x.condcode == 'FIX':
        return 'FIX'
    return x.condcode + str(x.itemnum)


def get_nelson_scores(t, pending=0, processed=0, closed=0):
    assert len(t.ch) < 3, 'Non-binary tree. %s' % t
    n_pending = []
    n_closed = []
    opennodes = []
    if len(t.ch):
        for i, ch in enumerate(t.ch):
            _n_pending, _n_closed, _opennodes = get_nelson_scores(ch, pending=pending, processed=processed, closed=closed + ((i==1) and (len(t.ch[0].ch) > 1)))
            pending += len(_n_pending)
            if len(_n_pending) > 1:
                processed += len(_n_pending)
            n_pending += _n_pending
            n_closed += _n_closed
            opennodes += _opennodes
    else:
        n_pending = [pending - processed + 1]
        n_closed = [closed]
        opennodes = [pending - processed + 1 + closed]

    return n_pending, n_closed, opennodes


# Collect item-level metrics

ling_preds = [
    'word',
    'wlen',
    'docid',
    'sentid',
    'sentpos',
    'cond',
    'condcode',
    'conlen',
    'dlt',
    'dltc',
    'dltv',
    'dltm',
    'dltcv',
    'dltcm',
    'dltvm',
    'dltcvm',
    'dlts',
    'unigramsurp',
    'fwprob5surp',
    'totsurp',
    'startembdAny',
    'endembdAny',
    'embddepthAny',
    'embdlen',
    'noF',
    'noFlen',
    'noFlenlog1p',
    'noFdr',
    'noFdrv',
    'yesJ',
    'embddepthMin',
    'opennodes',
    'nmerged'
]

ling_preds_full = ling_preds + ['PMI']

ling_preds_nojab = [
    'dlt',
    'dltc',
    'dltv',
    'dltm',
    'dltcv',
    'dltcm',
    'dltvm',
    'dltcvm',
    'dlts',
    'unigramsurp',
    'fwprob5surp',
    'totsurp',
    'startembdAny',
    'endembdAny',
    'embddepthAny',
    'embdlen',
    'noF',
    'noFlen',
    'noFlenlog1p',
    'noFdr',
    'noFdrv',
    'yesJ',
    'embddepthMin',
    'opennodes',
    'nmerged'
]

ling_baselines = [
    ['opennodes'],
    ['nmerged'],
    ['dlts'],
    ['dltcvm'],
    ['fwprob5surp'],
    ['totsurp']
]
ling_baselines.append([x[0] for x in ling_baselines])

itemmeasures = pd.read_csv(
    'ling_preds/conlen2fmri.wsj02to21-gcg15-nol-prtrm-3sm-synproc-+c_+u_+b5000_parsed.dlt.lc.unigram.5-kenlm.all-itemmeasures',
    sep=' '
)
itemmeasures['noFlenlog1p'] = np.log(itemmeasures['noFlen'].values + 1)

pending = []
closed = []
opennodes = []
nmerged = []
t = Tree()
with open('ling_preds/conlenc.gold.linetrees', 'r') as f:
    for line in f.readlines():
        t.read(line)
        t.collapseUnary()
        _pending, _closed, _opennodes = get_nelson_scores(t)
        opennodes += _opennodes
        _pending.append(0)
        _nmerged = [max(x - y + 1, 0) for x, y in zip(_pending[:-1], _pending[1:])]
        nmerged += _nmerged

opennodes_df = np.ones(len(itemmeasures), dtype=int)
opennodes_df[:len(opennodes)] = opennodes
itemmeasures['opennodes'] = opennodes_df
nmerged_df = np.zeros(len(itemmeasures), dtype=int)
nmerged_df[:len(opennodes)] = nmerged
itemmeasures['nmerged'] = nmerged_df
itemmeasures = itemmeasures[ling_preds]
itemmeasures.loc[itemmeasures.cond == 'JAB', ling_preds_nojab] = 0
itemmeasures = pd.merge(itemmeasures, df_pmi[['PMI', 'docid', 'sentpos', 'cond', 'conlen', 'condcode']],
                        how='left', on=['docid', 'sentpos', 'cond', 'conlen', 'condcode'])
itemmeasures = itemmeasures.fillna(0.)
itemmeasures['itempos'] = itemmeasures.groupby('docid').cumcount() + 1
itemmeasures['chunkpos'] = (itemmeasures.groupby('docid').cumcount()) % itemmeasures['conlen'] + 1
itemmeasures['chunkstart'] = ((itemmeasures['chunkpos'] - itemmeasures['chunkpos'].shift()) != 1)
itemmeasures['chunkid'] = itemmeasures['chunkstart'].cumsum()

itemmeans = itemmeasures[ling_preds_nojab + ['docid', 'cond', 'condcode', 'conlen']] \
    .groupby(['docid', 'cond', 'condcode', 'conlen']).mean().reset_index()


# Plot overall statistics

FONT_SIZE = 16
TICK_SIZE = 14
bar_width = 0.8
capthick=1
capsize=2
x_sep = 1
x_width = 1
n_network = 1
legend_size = 10
rotation = 30

colors_bycond = [
    # C
    (240, 202, 0, 255),
    (246, 158, 0, 255),
    (250, 122, 0, 255),
    (252, 90, 0, 255),
    (252, 66, 0, 255),
    (253, 42, 0, 255),

    # NC
    (39, 196, 246, 255),
    (7, 124, 206, 255),

    # JAB
    (222, 153, 255, 255),
    (175, 133, 238, 255),
    (160, 82, 202, 255),
]
colors_bycond = [tuple([float(x) / 255 for x in y]) for y in colors_bycond]
colors_bylen = [
    # C
    (252, 66, 0, 255),

    # NC
    (7, 124, 206, 255),

    # JAB
    (160, 82, 202, 255),
]
colors_bylen = [tuple([float(x) / 255 for x in y]) for y in colors_bylen]
orangered = colors.to_rgba('orangered')
dodgerblue = colors.to_rgba('dodgerblue')
gray = colors.to_rgba('gray')
gray = tuple([x if i < 3 else 0.4 for i, x in enumerate(gray)])

plt.rcParams.update({'font.size': FONT_SIZE, 'xtick.labelsize': TICK_SIZE, 'ytick.labelsize': TICK_SIZE})
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"

pos = pd.read_csv('ling_preds/pos.csv')
pos_means = {}
lb = {}
ub = {}
for pos_tag, v in pos.groupby('Part of Speech'):
    v = v.sort_values('Length')
    pos_means[pos_tag] = v['Mean'].values
    lb[pos_tag] = v['Mean'].values - v['2.5%'].values
    ub[pos_tag] = v['97.5%'].values - v['Mean'].values

fig = plt.figure(figsize=(10,3.1))
h = [Size.Fixed(1.0), Size.Fixed(6)]
v = [Size.Fixed(0.6), Size.Fixed(2.5)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax = fig.add_axes(
    divider.get_position(),
    axes_locator=divider.new_locator(nx=1, ny=1)
)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.tick_params(labelleft='on', labelbottom='on')
ax.yaxis.set_ticks_position('none')
ax.xaxis.set_ticks_position('bottom')

hatches = [
    '///',
    '...',
    '++',
    'oo'
]

left = np.zeros(6)
for i, pos_tag in enumerate(['Adjective/Adverb', 'Verb', 'Noun', 'Function Word']):
    ax.barh(
        np.arange(5,-1,-1),
        pos_means[pos_tag],
        left=left,
        # color=[get_color(x, c) for x in names],
        color='white',
        edgecolor='gray',
        hatch=hatches[i],
        lw=1.5,
        label=pos_tag
    )
    left += pos_means[pos_tag]

ax.set_yticks(np.arange(5,-1,-1))
ax.set_yticklabels(['c01', 'c02', 'c03', 'c04', 'c06', 'c12'])
ax.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')
ax.set_xlabel('Proportion of words')

if not os.path.exists('plots'):
    os.makedirs('plots')

plt.savefig('plots/pos_distribution.png')

fig = plt.figure(figsize=((3, 9./4)))
h = [Size.Fixed(0.5), Size.Fixed(2.5)]
v = [Size.Fixed(0.5), Size.Fixed(7./4)]
divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax = fig.add_axes(
    divider.get_position(),
    axes_locator=divider.new_locator(nx=1, ny=1)
)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.tick_params(labelleft='on', labelbottom='on')
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

ax.barh(
    np.arange(5,-1,-1),
    np.arange(6),
    color=colors_bycond,
    edgecolor='none',
    lw=1.5,
    label=['c12', 'c06', 'c04', 'c03', 'c02', 'c01'],
)

ax.set_yticks(np.arange(5,-1,-1) * x_width)
ax.set_yticklabels(['c01', 'c02', 'c03', 'c04', 'c06', 'c12'], rotation=rotation, ha='right')

if not os.path.exists('plots'):
    os.makedirs('plots')

plt.savefig('plots/items.pdd.png')
plt.close('all')

for ling_pred in ling_baselines:
    if len(ling_pred) == 1:
        ling_pred = ling_pred[0]
        df = itemmeans[itemmeans.cond == 'C']
        ling_mean = []
        ling_err = []
        for x in [1, 2, 3, 4, 6, 12]:
            ling_mean.append(df[df.conlen == x][ling_pred].mean())
            ling_err.append(df[df.conlen == x][ling_pred].sem())

        fig = plt.figure(figsize=((3, 9./4)))
        h = [Size.Fixed(0.5), Size.Fixed(2.5)]
        v = [Size.Fixed(0.5), Size.Fixed(7./4)]
        divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
        ax = fig.add_axes(
            divider.get_position(),
            axes_locator=divider.new_locator(nx=1, ny=1)
        )
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.tick_params(labelleft='on', labelbottom='on')
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        ax.barh(
            np.arange(5,-1,-1),
            ling_mean,
            color=colors_bycond,
            edgecolor='none',
            lw=1.5,
            label=['c12', 'c06', 'c04', 'c03', 'c02', 'c01'],
        )

        for j, x in enumerate(['c12', 'c06', 'c04', 'c03', 'c02', 'c01']):
            ax.errorbar(
                ling_mean[j:j + 1],
                [5-j],
                xerr=ling_err[j:j + 1],
                fmt='none',
                ecolor='black',
                lw=2,
                capthick=capthick,
                capsize=capsize
            )

        ax.set_yticks(np.arange(5,-1,-1) * x_width)
        ax.set_yticklabels(['c01', 'c02', 'c03', 'c04', 'c06', 'c12'], rotation=rotation, ha='right')

        plt.savefig('plots/items.%s.png' % ling_pred)
        plt.close('all')
