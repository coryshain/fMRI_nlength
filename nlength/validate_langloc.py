import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import fdrcorrection

try:
    with open('data_path.txt', 'r') as f:
        base_path = os.path.join(f.read().strip(), 'main')
except FileNotFoundError:
    sys.stderr.write('Data path not set. Run `python -m nlength.set_data_path` before running any other scripts.\n')
    sys.stderr.flush()
    exit()

expts = [
    '6words_n20',
    'new_subjects_n15',
    'nlength_con_n16',
    'old_subjects_n25',
]

contrasts = {}
for expt in expts:
    for localizer in ('lang', 'alice'):
        path = '%s%s/func_parcels/mROI_%sEFFECT_%sLOC/spm_ss_mROI_data.csv' % (base_path, expt, localizer, localizer)
        if os.path.exists(path):
            df = pd.read_csv(path)
            df = df.drop_duplicates(['Subject', 'ROI', 'Effect'])
            for ROI in range(1, 7):
                if localizer == 'lang':
                    df_a = df[(df.ROI == ROI) & (df.Effect == 'S')].sort_values('Subject')
                    df_b = df[(df.ROI == ROI) & (df.Effect == 'N')].sort_values('Subject')
                else:
                    df_a = df[(df.ROI == ROI) & (df.Effect == 'I')].sort_values('Subject')
                    df_b = df[(df.ROI == ROI) & (df.Effect == 'D')].sort_values('Subject')
                contrast = df_a.EffectSize.values - df_b.EffectSize.values
                if not ROI in contrasts:
                    contrasts[ROI] = []
                contrasts[ROI].append(contrast)

for ROI in contrasts:
    contrasts[ROI] = np.concatenate(contrasts[ROI])

t = []
p = []
df = []
d = []
for ROI in contrasts:
    x = contrasts[ROI]
    out = ttest_1samp(x, 0)
    _t = out.statistic
    _p = out.pvalue
    _df = out._df
    _d = x.mean() / x.std()

    t.append(_t)
    p.append(_p)
    df.append(_df)
    d.append(_d)

p = fdrcorrection(p, method='negcorr')[1]

for ROI, _t, _p, _df, _d in zip(contrasts.keys(), t, p, df, d):
    print('ROI: %s | t: %s | p: %s | df: %s | d: %s' % (ROI, _t, _p, _df, _d))


