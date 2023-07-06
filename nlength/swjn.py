import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import sem, ttest_1samp
from statsmodels.stats.multitest import fdrcorrection
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size

from nlength.signif_table import correct_p, get_network_fdr, get_stars

def map_effect(x):
    y = x.split('_')[-1]
    if y == '12c':
        return 'S'
    if y == '1c':
        return 'W'
    if y == 'jab12c':
        return 'J'
    if y == 'jab1c':
        return 'N'
    return x

try:
    with open('data_path.txt', 'r') as f:
        base_path = f.read().strip()
except FileNotFoundError:
    sys.stderr.write('Data path not set. Run `python -m nlength.set_data_path` before running any other scripts.\n')
    sys.stderr.flush()
    exit()

df_fed_exp1 = pd.read_csv(os.path.join(base_path, 'fed10_data/SWJNV1_results.csv'), sep=',\s?').sort_values(['Subject', 'ROI'])
df_fed_exp2 = pd.read_csv(os.path.join(base_path, 'fed10_data/SWJNV2_results.csv'), sep=',\s?').sort_values(['Subject', 'ROI'])
df_curr = []
for subj_set in ('old_subjects_n25', 'new_subjects_n15'):
    df_curr.append(pd.read_csv(os.path.join(base_path, 'main', subj_set, 'func_parcels', 'mROI_NlengthEFFECT_langLOC',
                     'spm_ss_mROI_data.details.EffectSize.csv')))
    if subj_set == 'old_subjects_n25':
        df_curr.append(pd.read_csv(os.path.join(base_path, 'main', subj_set, 'func_parcels', 'mROI_NlengthEFFECT_langrun1LOC',
                     'spm_ss_mROI_data.details.EffectSize.csv')))
df_curr = pd.concat(df_curr, axis=0)
df_curr = df_curr[~df_curr.Effect.str.contains('-')]
df_curr.Effect = df_curr.Effect.apply(map_effect)

fROIs = {
    1: 'LIFGorb',
    2: 'LIFG',
    3: 'LMFG',
    4: 'LAntTemp',
    5: 'LPostTemp',
    6: 'LAngG'
}


# Report tests of S > N

df = []
for i in range(1, 7):
    for exp, _df in zip(('fed_exp1', 'fed_exp2', 'current'), (df_fed_exp1, df_fed_exp2, df_curr)):
        contrasts = {}
        contrasts['S'] = _df[(_df.ROI == i) & (_df.Effect == 'S')].EffectSize.values
        contrasts['W'] = _df[(_df.ROI == i) & (_df.Effect == 'W')].EffectSize.values
        contrasts['J'] = _df[(_df.ROI == i) & (_df.Effect == 'J')].EffectSize.values
        contrasts['N'] = _df[(_df.ROI == i) & (_df.Effect == 'N')].EffectSize.values
        contrasts['S_v_W'] = contrasts['S'] - contrasts['W']
        contrasts['S_v_N'] = contrasts['S'] - contrasts['N']
        contrasts['J_v_N'] = contrasts['J'] - contrasts['N']
        contrasts['S_v_W_v_J_v_N'] = contrasts['S_v_W'] - contrasts['J_v_N']
        for contrast in contrasts:
            vals = contrasts[contrast]
            t, p = ttest_1samp(vals, 0.)
            beta = vals.mean()
            se = sem(vals)
            d = beta / vals.std(ddof=1)
            df.append((exp, contrast, fROIs[i], beta, se, t, p, d))

df = pd.DataFrame(df, columns=['experiment', 'contrast', 'fROI', 'beta', 'se', 't', 'p', 'd'])
df.fROI = pd.Categorical(
    df.fROI,
    [
        'all',
        'LIFGorb',
        'LIFGtri',
        'LIFG',
        'LMFG',
        'LTP',
        'LaSTS',
        'LAntTemp',
        'LpSTS',
        'LPostTemp',
        'LTPJ',
        'LAngG',
        'RIFGorb',
        'RIFG',
        'RMFG',
        'RAntTemp',
        'RPostTemp',
        'RAngG',
        'IFGorb',
        'IFG',
        'MFG',
        'AntTemp',
        'PostTemp',
        'AngG',
        'LIFGorb_v_LAngG',
        'LIFG_v_LAngG',
        'LMFG_v_LAngG',
        'LAntTemp_v_LAngG',
        'LPostTemp_v_LAngG',
        'LIFGorb_v_LPostTemp',
        'LIFG_v_LPostTemp'
    ]
)
df['network_fdr'] = df[['fROI']].apply(get_network_fdr, axis=1)
out_p = df.groupby(['experiment', 'network_fdr', 'contrast']) \
        .apply(correct_p)[['experiment', 'network_fdr', 'contrast', 'fROI', 'p_fdr']]
df = pd.merge(df, out_p, on=['experiment', 'network_fdr', 'contrast', 'fROI'])
df['signif'] = df['p_fdr'].apply(get_stars)
df.beta = df.beta.astype(float).round(2)
df.se = df.se.astype(float).round(2)
df.t = df.t.astype(float).round(2)
df.p = np.maximum(df.p.astype(float), 0.001).round(3)
df.p_fdr = np.maximum(df.p_fdr.astype(float), 0.001).astype(float).round(3)
df = df.sort_values(['experiment', 'contrast', 'fROI'])
df.to_csv('swjn.csv')


# Plot

font_size = 16
tick_size = 14
legend_size = 10
plt.rcParams.update({'font.size': font_size, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size})
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"

fig = plt.figure(figsize=((13, 6.2)))
h = [Size.Fixed(0.7), Size.Fixed(10.3)]
v = [Size.Fixed(1.5), Size.Fixed(4.1)]
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
ax.axhline(y=0, lw=1, c='gray', alpha=1, zorder=2)

color = [
    # S
    (137, 69, 246, 255),

    # W
    (14, 60, 245, 255),

    # J
    (235, 63, 37, 255),

    # N
    (192, 192, 192, 255)
]
color = [tuple([float(x) / 255 for x in y]) for y in color]

bar_width = 1. / 12 * 0.7
contrasts = ['S', 'W', 'J', 'N']

for i in range(1, 6):
    ax.axvline(x=i - (1.5 * bar_width), lw=1, c='gray', alpha=1, zorder=2)

r_base = np.arange(6)
conds = ('S', 'W', 'J', 'N')
dfs = (df_fed_exp1, df_fed_exp2, df_curr)
for i in range(12):
    j = i % 4
    k = i // 4
    cond = conds[j]
    _df = dfs[k]
    estimates = []
    errors = []
    estimates = _df[(_df.Effect == cond) & (_df.ROI <= 6)].groupby('ROI')['EffectSize'].mean()
    errors = _df[(_df.Effect == cond) & (_df.ROI <= 6)].groupby('ROI')['EffectSize'].sem()

    r = r_base + i * bar_width + (i // 4) * 0.1

    ax.bar(
        r,
        estimates,
        color=color[j],
        width=bar_width,
        label=cond if k == 0 else None,
        linewidth=2,
        # linestyle='solid' if (i % 2 == 0) else 'dashed'
    )

    ax.errorbar(
        r,
        estimates,
        yerr=errors,
        fmt='none',
        ecolor=color[j],
        capsize=2,
        capthick=2,
        linewidth=2
    )

ax.legend(loc='center left', ncol=1, bbox_to_anchor=(1, 0.5))

ax.set_xticks(np.arange(0, 5.9, 0.3333333) + bar_width * 1.5)
ax.set_xticklabels(['F. et al (2010) Exp1', 'F. et al (2010) Exp2', 'Current Study'] * 6, rotation=45, ha='right')
# ax.set_ylabel('BOLD')
ax.set_ylim(-1.37, 6)

if not os.path.exists('plots'):
    os.makedirs('plots')

plt.savefig('plots/SWJN.png')



