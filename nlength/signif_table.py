import sys
import os
import re
import argparse
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection


def correct_p(x):
    p = fdrcorrection(x['p'].values, method='negcorr')[1]
    x['p_fdr'] = p
    return x


def get_stars(x):
    if x > 0.1:
        return ''
    if x > 0.05:
        return '.'
    if x > 0.01:
        return '*'
    if x > 0.001:
        return '**'
    return '***'


def get_network_fdr(x):
    if x.fROI != 'all':
        return 'Ind'
    return x.fROI


# Thanks to Daniel Sparks on StackOverflow for this one (post available at
# http://stackoverflow.com/questions/5084743/how-to-print-pretty-string-output-in-python)
def pretty_table(row_collection, key_list, field_sep=' '):
  return '\n'.join([field_sep.join([str(row[col]).ljust(width)
    for (col, width) in zip(key_list, [max(map(len, column_vector))
      for column_vector in [ [v[k]
        for v in row_collection if k in v]
          for k in key_list ]])])
            for row in row_collection])


def csv_table(row_collection, key_list):
    return '\n'.join([','.join(z) for z in ([[str(x[y]) for y in key_list] for x in row_collection])])


def compute_row(path, stars=True):
    experiment = os.path.basename(path).split('.')[0].split('_')[1]
    contrast, fROI = path.split('.')[-5:-3]
    conv = 't'
    sing = 'f'
    beta = '---'
    se = '---'
    t = '---'
    p = '---'
    betweengroups = 'betweengroups' in path

    with open(path, 'r') as f:
        in_lrt = False
        in_full = False
        in_fixef = False
        for line in f:
            if 'isSingular' in line:
                sing = 't'
            if 'failed to converge' in line:
                conv = 'f'
            if not line.strip():
                if in_full and in_fixef:
                    in_full = False
                    in_fixef = False
            if line.startswith('Full model'):
                in_full = True
            elif line.startswith('Fixed effects:') or line.startswith('Coefficients:'):
                in_fixef = True
            elif betweengroups and line.startswith('experiment') or (not betweengroups and line.startswith('(Intercept)')):
                if in_full and in_fixef:
                    _, beta, se, t = line.strip().split()[:4]
            elif line.startswith('isC'):
                if in_full and in_fixef:
                    _, beta, se, t = line.strip().split()[:4]
            elif line.strip().startswith('npar') or line.strip().startswith('Res.Df'):
                in_lrt = True
            elif line.startswith('m_full') or line.startswith('2'):
                if in_lrt:
                    line_parts = line.strip().split()
                    if '<' in line_parts:
                        line_parts.remove('<')
                    if len(line_parts) > 8:
                        p = float(line_parts[8])
                    elif len(line_parts) > 6:
                        p = float(line_parts[6])
                    else:
                        p = 1
                    if stars:
                        stars = 0
                        if p < 0.05:
                            stars += 1
                        if p < 0.01:
                            stars += 1
                        if p < 0.001:
                            stars += 1
                    p = '%.3e' % p
                    if stars:
                        p += ('*' * stars)
                    in_lrt = False

    return {
        'experiment': experiment,
        'contrast': contrast,
        'fROI': fROI,
        'conv': conv,
        'sing': sing,
        'beta': beta,
        'se': se,
        't': t,
        'p': p
    }


if __name__ == '__main__':
    argparser = argparse.ArgumentParser('''
    Get table of significance values for conlen tests.
    ''')
    args = argparser.parse_args()
    cols = ['parcel_set', 'experiment', 'contrast', 'fROI', 'conv', 'sing', 'beta', 'se', 't', 'p']
    rows = []
    for parcel_set in ['evlab', 'PDD', 'RH', 'betweenhemispheres', 'PDDanat']:
        for experiment in ['nlength1', 'nlength2', '6words', 'betweengroups']:
            directory = 'lrt'
            paths = sorted([os.path.join(directory, x) for x in os.listdir(directory)
                    if (x.startswith('%s_%s' % (parcel_set, experiment)) and x.endswith('.lrt.summary.txt'))])
            for path in paths:
                row = compute_row(path, stars=False)
                row['parcel_set'] = parcel_set
                rows.append(row)

    rows = sorted(rows, key=lambda x: (x['parcel_set'], x['experiment'], len(x['contrast']), x['contrast'].split('!')[0]))

    df = pd.DataFrame(rows, columns=cols)
    df['p'] = df['p'].astype(float)
    df['network_fdr'] = df[['fROI']].apply(get_network_fdr, axis=1)
    df_p = df.groupby(['parcel_set', 'experiment', 'network_fdr', 'contrast']) \
        .apply(correct_p)[['parcel_set', 'experiment', 'network_fdr', 'contrast', 'fROI', 'p_fdr']].reset_index(drop=True)
    df = pd.merge(df, df_p, on=['parcel_set', 'experiment', 'network_fdr', 'contrast', 'fROI'])
    df['signif'] = df['p_fdr'].apply(get_stars)
    df.beta = df.beta.astype(float).round(2)
    df.se = df.se.astype(float).round(2)
    df.t = df.t.astype(float).round(2)
    df.p = np.maximum(df.p.astype(float), 0.001).round(3)
    df.p_fdr = np.maximum(df.p_fdr.astype(float), 0.001).astype(float).round(3)
    df.parcel_set = pd.Categorical(df.parcel_set, ['evlab', 'PDD', 'RH', 'betweenhemispheres', 'PDDanat'])
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
    df = df.sort_values(['parcel_set', 'experiment', 'contrast', 'fROI'])
    df.to_csv('signif.csv', index=False)
