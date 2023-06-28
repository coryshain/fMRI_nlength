import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


evlab_parcels =  [
    'LIFGorb',
    'LIFG',
    'LMFG',
    'LAntTemp',
    'LPostTemp',
    'LAngG',
    'RIFGorb',
    'RIFG',
    'RMFG',
    'RAntTemp',
    'RPostTemp',
    'RAngG'
]

PDD_parcels = [
    'LIFGorb',
    'LIFGtri',
    'LTP',
    'LaSTS',
    'LpSTS',
    'LTPJ',
]

fROIs = {
    'evlab': evlab_parcels,
    'PDD': PDD_parcels,
    'RH': evlab_parcels,
}

networks = {
    'evlab': [
        'LIFGorb',
        'LIFG',
        'LMFG',
        'LAntTemp',
        'LPostTemp',
        'LAngG',
    ],
    'PDD': [
        'LIFGorb',
        'LIFGtri',
        'LTP',
        'LaSTS',
        'LpSTS',
        'LTPJ',
    ],
    'RH': [
        'RIFGorb',
        'RIFG',
        'RMFG',
        'RAntTemp',
        'RPostTemp',
        'RAngG',
    ],
}

parcel_set_map = {
    'evlab': 'func_parcels',
    'PDD': 'PDD_parcels',
    'RH': 'func_parcels',
}
experiments = (1, 2)
base_path = '../../data/fMRI_nlength/casto'

LENGTH2X = {
    1: 0.,
    2: 1.,
    3: 2.,
    4: 3.,
    5: 3.5,
    6: 4.,
    8: 4.33,
    10: 4.66,
    12: 5.
}

def length2x(x):
    x = np.array(x)
    f = np.vectorize(lambda x: LENGTH2X[int(x)])
    return f(x)

for parcel_set in parcel_set_map:
    parcel_set_path = parcel_set_map[parcel_set]
    for experiment in experiments:
        if experiment == 1:
            df = pd.read_csv(os.path.join(base_path, 'no_npmod_FINAL', 'nlength_con_n16', parcel_set_path, 'mROI_NlengthEFFECT_langLOC',
                             'spm_ss_mROI_data.details.EffectSize.csv'))
            df = df[df.Subject != '430_FED_20170523b_3T2_PL2017']  # Drop repeated session by subject 430
            lengths = [1, 2, 4, 6, 12]  # Length 3 condition missing
        else:
            df = []
            for subj_set in ('old_subjects_n25', 'new_subjects_n15'):
                df.append(pd.read_csv(os.path.join(base_path, 'no_npmod_FINAL', subj_set, parcel_set_path, 'mROI_NlengthEFFECT_langLOC',
                                 'spm_ss_mROI_data.details.EffectSize.csv')))
                if subj_set == 'old_subjects_n25':
                    df.append(pd.read_csv(os.path.join(base_path, 'no_npmod_FINAL', subj_set, parcel_set_path, 'mROI_NlengthEFFECT_langrun1LOC',
                                 'spm_ss_mROI_data.details.EffectSize.csv')))
            df = pd.concat(df, axis=0)
            lengths = [1, 2, 3, 4, 6, 12]

        plot_basis = length2x(lengths)

        xtick_pos = plot_basis
        xtick_labels = [str(x) for x in lengths]

        plot_path = 'plots'
        contrast_path = 'contrasts'

        if not os.path.exists(plot_path):
            os.makedirs(plot_path)

        if not os.path.exists(contrast_path):
            os.makedirs(contrast_path)

        df.ROI = df.ROI.apply(lambda x: fROIs[parcel_set][x-1])
        df = df[df.ROI.isin(networks[parcel_set])]
        df = df[~df.Effect.str.contains('-')]
        df['StimType'] = np.zeros_like(df.Effect)
        df.StimType[df.Effect.str.contains('jab')] = 'J'
        df.StimType[df.Effect.str.contains('nc')] = 'N'
        df.StimType[(~df.Effect.str.contains('nc')) & (~df.Effect.str.contains('jab'))] = 'C'
        df['nlength'] = df.Effect.str.extract('(\d+)').astype(int)

        out = []


        for ROI in networks[parcel_set]:
            # Plot
            plt.gca().spines['top'].set_visible(False)
            plt.gca().spines['right'].set_visible(False)
            plt.gca().spines['bottom'].set_visible(False)
            plt.gca().spines['left'].set_visible(True)
            plt.gca().tick_params(labelleft='on', labelbottom='on')
            plt.gca().yaxis.set_ticks_position('left')
            plt.gca().xaxis.set_ticks_position('none')
            # ax.grid(b=True, which='major', axis='y', ls='--', lw=.5, c='k', alpha=.3)
            plt.gca().axhline(y=0, lw=1, c='gray', alpha=1)

            _df = df[df.StimType == 'C']
            clens = [1, 2, 3, 4, 6, 12]
            means = []
            errs = []
            D_C = []
            subjects = None
            for i, clen in enumerate(clens):
                d = _df[(_df.nlength == clen) & (_df.ROI == ROI)]
                if len(d.values):
                    d = d.sort_values('Subject')
                    if subjects is None:
                        subjects = d.Subject.values
                    d = d.EffectSize
                    m = d.mean()
                    means.append(m)
                    sem = d.sem()
                    errs.append(sem)
                    D_C.append(d.values)
            D_C = np.stack(D_C, axis=1)

            b = np.linalg.lstsq(np.stack([np.ones_like(means), plot_basis], axis=1), D_C.T)[0]
            NLenC = b[1]

            xline = np.linspace(plot_basis.min(), plot_basis.max(), 500)
            X = np.stack([np.ones_like(xline), xline], axis=1)
            yline = np.dot(X, b).mean(axis=-1)

            plt.errorbar(
                plot_basis,
                means,
                yerr=errs,
                fmt='ro',
                linestyle='none',
                ecolor='red',
                lw=2,
                capsize=0,
                label='normal'
            )
            plt.plot(
                xline,
                yline,
                linestyle='dashed',
                color='red',
            )

            if experiment == 2:
                _df = df[df.StimType == 'J']
                clens = [1, 4, 12]
                means = []
                errs = []
                D_J = []
                for i, clen in enumerate(clens):
                    d = _df[(_df.nlength == clen) & (_df.ROI == ROI)]
                    d = d.sort_values('Subject')
                    d = d.EffectSize
                    m = d.mean()
                    means.append(m)
                    sem = d.sem()
                    errs.append(sem)
                    D_J.append(d.values)
                D_J = np.stack(D_J, axis=1)

                b = np.linalg.lstsq(np.stack([np.ones_like(means), length2x(clens)], axis=1), D_J.T)[0]
                NLenJ = b[1]
                xline = np.linspace(plot_basis.min(), plot_basis.max(), 500)
                X = np.stack([np.ones_like(xline), xline], axis=1)
                yline = np.dot(X, b).mean(axis=-1)

                plt.errorbar(
                    [0, 3, 5],
                    means,
                    yerr=errs,
                    fmt='bs',
                    linestyle='none',
                    ecolor='blue',
                    lw=2,
                    capsize=0,
                    label='normal'
                )
                plt.plot(
                    xline,
                    yline,
                    linestyle='dashed',
                    color='blue'
                )

                _df = df[df.StimType == 'N']
                clens = [3, 4]
                means = []
                errs = []
                D_N = []
                for i, clen in enumerate(clens):
                    d = _df[(_df.nlength == clen) & (_df.ROI == ROI)]
                    d = d.sort_values('Subject')
                    d = d.EffectSize
                    m = d.mean()
                    means.append(m)
                    sem = d.sem()
                    errs.append(sem)
                    D_N.append(d.values)
                D_N = np.stack(D_N, axis=1)

                b = np.linalg.lstsq(np.stack([np.ones_like(means), length2x(clens)], axis=1), D_N.T)[0]
                NLenN = b[1]

                xline = np.linspace(plot_basis.min(), plot_basis.max(), 500)
                X = np.stack([np.ones_like(xline), xline], axis=1)
                yline = np.dot(X, b).mean(axis=-1)

                plt.errorbar(
                    [2, 3],
                    means,
                    yerr=errs,
                    fmt='mx',
                    linestyle='none',
                    ecolor='m',
                    lw=2,
                    capsize=0,
                    label='normal'
                )
                plt.plot(
                    xline,
                    yline,
                    linestyle='dashed',
                    color='m'
                )

            else:
                D_J = NLenJ = D_N = NLenN = None

            plt.subplots_adjust(left=0.3)
            plt.xlim(plot_basis.min() - 0.2, plot_basis.max() + 0.2)
            plt.ylim((-0.68, 3.5))

            plt.xticks(xtick_pos, labels=xtick_labels)

            plt.gcf().set_size_inches((2, 3))
            plt.savefig(os.path.join('plots', '%s_nlength%s_%s_plot.png' % (parcel_set, experiment, ROI)), dpi=300)
            plt.close('all')

            # Contrasts
            columns = ['C%02d' % x for x in lengths]
            if experiment == 2:
                columns += ['J%02d' % x for x in [1, 4, 12]] + ['N%02d' % x for x in [3, 4]]

            _out = [D_C]
            if experiment == 2:
                _out += [D_J, D_N]
            _out = pd.DataFrame(np.concatenate(_out, axis=1), columns=columns)
            _out['Subject'] = subjects
            _out['fROI'] = ROI

            C = np.dot(D_C, np.ones_like(plot_basis) / 6.)
            if experiment == 2:
                C1412 = np.dot(D_C, [0.333, 0, 0, 0.333, 0, 0.333])
                C126 = np.dot(D_C, [0.33, 0.33, 0, 0, 0.33, 0])
                C34 = np.dot(D_C, [0, 0, 0.5, 0.5, 0, 0])
                J = np.dot(D_J, [0.333, 0.333, 0.333])
                N = np.dot(D_N, [0.5, 0.5])
            else:
                C1412 = np.dot(D_C, [0.333, 0, 0.333, 0, 0.333])
                C126 = np.dot(D_C, [0.33, 0.33, 0, 0.33, 0])
                C34 = J = N = None

            S_v_W = D_C[:,-1] - D_C[:,0]
            if experiment == 2:
                S_v_N = D_C[:,-1] - D_J[:,0]
                J_v_W = D_J[:,-1] - D_C[:,0]
                J_v_N = D_J[:,-1] - D_J[:,0]
            else:
                S_v_N = J_v_W = J_v_N = None

            _out['C'] = C
            _out['C1412'] = C1412
            _out['C126'] = C126
            if experiment == 2:
                _out['C34'] = C34
                _out['N'] = N
                _out['J'] = J
                _out['C_v_J'] = C1412 - J
                _out['C_v_N'] = C34 - N

            _out['NLenC'] = NLenC
            if experiment == 2:
                _out['NLenJ'] = NLenJ
                _out['NLenN'] = NLenN

                _out['NLenC_v_NLenJ'] = NLenC - NLenJ
                _out['NLenC_v_NLenN'] = NLenC - NLenN
                _out['NLenJ_v_NLenN'] = NLenJ - NLenN

            _out['S_v_W'] = S_v_W
            if experiment == 2:
                _out['S_v_N'] = S_v_N
                _out['J_v_W'] = J_v_W
                _out['J_v_N'] = J_v_N

            out.append(_out)

        out = pd.concat(out, axis=0)
        out.to_csv(os.path.join('contrasts', '%s_nlength%s_contrasts.csv' % (parcel_set, experiment)), index=False)


    # 6 words experiment

    paths = [
        os.path.join(base_path, 'no_npmod_FINAL', '6words_n20', parcel_set_path, 'mROI_6wordsEFFECT_langLOC',
                     'spm_ss_mROI_data.details.EffectSize.csv'),
        os.path.join(base_path, 'no_npmod_FINAL', '6words_n20', parcel_set_path, 'mROI_6wordsEFFECT_aliceLOC',
                     'spm_ss_mROI_data.details.EffectSize.csv'),
    ]
    lengths = [1, 2, 3, 4, 5, 6, 8, 10, 12]
    wl24_lengths = [1, 2, 3, 4, 6, 8, 12]
    wl24_plot_basis = length2x(wl24_lengths)
    wl30_lengths = [1, 2, 3, 5, 6, 10]
    wl30_plot_basis = length2x(wl30_lengths)

    xtick_pos = length2x(lengths)
    xtick_labels = [str(x) for x in lengths]

    df = [pd.read_csv(path) for path in paths]
    df = pd.concat(df, axis=0)

    plot_path = 'plots'
    contrast_path = 'contrasts'

    df.ROI = df.ROI.apply(lambda x: fROIs[parcel_set][x - 1])
    df = df[df.ROI.isin(networks[parcel_set])]
    df = df[~df.Effect.str.contains('-')]
    df['StimType'] = np.zeros_like(df.Effect)
    df.StimType[df.Effect.str.contains('jab')] = 'J'
    df.StimType[df.Effect.str.contains('nc')] = 'N'
    df.StimType[(~df.Effect.str.contains('nc')) & (~df.Effect.str.contains('jab'))] = 'C'
    df['nlength'] = df.Effect.str.extract('_?(\d+)[_cnj]').astype(int)
    df['experiment'] = df.Effect.apply(lambda x: '6words' if 'wl' in x else 'nlength')

    out = []

    for ROI in networks[parcel_set]:
        # Plot
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)
        plt.gca().spines['left'].set_visible(True)
        plt.gca().tick_params(labelleft='on', labelbottom='on')
        plt.gca().yaxis.set_ticks_position('left')
        plt.gca().xaxis.set_ticks_position('none')
        plt.gca().axhline(y=0, lw=1, c='gray', alpha=1)

        # 24 word list
        _df = df[df.Effect.str.endswith('24wl')]
        clens = wl24_lengths
        plot_basis = wl24_plot_basis
        means = []
        errs = []
        D_WL24 = []
        subjects = None
        for i, clen in enumerate(clens):
            d = _df[(_df.nlength == clen) & (_df.ROI == ROI)]
            if len(d.values):
                d = d.sort_values('Subject')
                if subjects is None:
                    subjects = d.Subject.values
                d = d.EffectSize
                m = d.mean()
                means.append(m)
                sem = d.sem()
                errs.append(sem)
                D_WL24.append(d.values)
        D_WL24 = np.stack(D_WL24, axis=1)

        b = np.linalg.lstsq(np.stack([np.ones_like(means), plot_basis], axis=1), D_WL24.T)[0]
        NLen24WL = b[1]

        xline = np.linspace(0, length2x(12), 500)
        X = np.stack([np.ones_like(xline), xline], axis=1)
        yline = np.dot(X, b).mean(axis=-1)

        plt.errorbar(
            plot_basis,
            means,
            yerr=errs,
            fmt='gv',
            linestyle='none',
            ecolor='green',
            lw=2,
            capsize=0,
            label='normal'
        )
        plt.plot(
            xline,
            yline,
            linestyle='dashed',
            color='green',
        )

        # 30 word list
        _df = df[df.Effect.str.endswith('30wl')]
        clens = wl30_lengths
        plot_basis = wl30_plot_basis
        means = []
        errs = []
        D_WL30 = []
        subjects = None
        for i, clen in enumerate(clens):
            d = _df[(_df.nlength == clen) & (_df.ROI == ROI)]
            if len(d.values):
                d = d.sort_values('Subject')
                if subjects is None:
                    subjects = d.Subject.values
                d = d.EffectSize
                m = d.mean()
                means.append(m)
                sem = d.sem()
                errs.append(sem)
                D_WL30.append(d.values)
        D_WL30 = np.stack(D_WL30, axis=1)

        b = np.linalg.lstsq(np.stack([np.ones_like(means), plot_basis], axis=1), D_WL30.T)[0]
        NLen30WL = b[1]

        xline = np.linspace(0, length2x(12), 500)
        X = np.stack([np.ones_like(xline), xline], axis=1)
        yline = np.dot(X, b).mean(axis=-1)

        plt.errorbar(
            plot_basis,
            means,
            yerr=errs,
            fmt='c^',
            linestyle='none',
            ecolor='c',
            lw=2,
            capsize=0,
            label='normal'
        )
        plt.plot(
            xline,
            yline,
            linestyle='dashed',
            color='c',
        )

        C_WL24 = np.dot(D_WL24, np.ones_like(wl24_lengths) / len(wl24_lengths))
        C_WL30 = np.dot(D_WL30, np.ones_like(wl30_lengths) / len(wl30_lengths))
        X = np.concatenate([D_WL24, D_WL30], axis=1)
        ncol = X.shape[1]
        C_6words = np.dot(X, np.ones(ncol) / ncol)
        C126_6words = np.dot(X, [0.166, 0.166, 0, 0, 0.166, 0, 0, 0.166, 0.166, 0, 0, 0.166, 0])
        steps = length2x(wl24_lengths + wl30_lengths)
        b = np.linalg.lstsq(np.stack([np.ones(ncol), steps], axis=1), X.T)[0]
        NLen6words = b[1]

        plt.subplots_adjust(left=0.3)
        plt.xticks(length2x([1, 2, 3, 4, 6, 12]), ['1', '2', '3', '4', '6', '12'])
        plt.xlim(-0.2, length2x(12) + 0.2)
        plt.ylim((-0.68, 3.5))

        plt.gcf().set_size_inches((2, 3))
        plt.savefig(os.path.join('plots', '%s_6words_%s_plot.png' % (parcel_set, ROI)), dpi=300)
        plt.close('all')

        columns = ['C%02d_WL24' % x for x in wl24_lengths]
        columns += ['C%02d_WL30' % x for x in wl30_lengths]

        _out = pd.DataFrame(np.concatenate([D_WL24, D_WL30], axis=1), columns=columns)
        _out['Subject'] = subjects
        _out['fROI'] = ROI


        _out['C_WL24'] = C_WL24
        _out['C_WL30'] = C_WL30
        _out['C_6words'] = C_6words
        _out['C126_6words'] = C126_6words
        _out['NLen24WL'] = NLen24WL
        _out['NLen30WL'] = NLen30WL
        _out['NLen6words'] = NLen6words

        out.append(_out)

    out = pd.concat(out, axis=0)
    out.to_csv(os.path.join('contrasts', '%s_6words_contrasts.csv' % parcel_set), index=False)


