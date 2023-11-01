#!/usr/bin/env Rscript

pr = function(x, buffer=NULL) {
    if (is.null(buffer)) {
        buffer = stderr()
    }
    cat(paste0(x, '\n'), file=buffer, append=TRUE)
}

library(lme4)


for (parcel_set in c('evlab', 'PDD', 'RH', 'PDDanat')) {
    for (experiment in c(1, 2, 3)) {
        if (experiment < 3) {
            df_name = paste0('nlength', experiment)
        } else {
            df_name = '6words'
        }
        df = read.table(paste0('contrasts/', parcel_set, '_', df_name, '_contrasts.csv'), header=TRUE, sep=',')
        prefix_out = paste0('glm/', parcel_set, '_', df_name)
        if (!dir.exists('glm')) {
            dir.create('glm', recursive = TRUE)
        }

        if (experiment == 1) {
            contrasts = c(
                'C',
                'NLenC',
                'S_v_W'
            )
        } else if (experiment == 2) {
            contrasts = c(
                'C',
                'J',
                'N',
                'C_v_J',
                'C_v_N',
                'NLenC',
                'NLenJ',
                'NLenN',
                'NLenC_v_NLenJ',
                'NLenC_v_NLenN',
                'NLenJ_v_NLenN',
                'S_v_W',
                'S_v_N',
                'J_v_W',
                'J_v_N',
                'S_v_W_v_J_v_N'
            )
        } else {
            contrasts = c(
                'C_WL24',
                'C_WL30',
                'C_6words',
                'NLen24WL',
                'NLen30WL',
                'NLen6words'
            )
        }

        for (contrast in contrasts) {
            for (fROI in c('all', unique(df$fROI))) {
                for (mtype in c('full', 'abl')) {
                    if (mtype == 'full') {
                        if (fROI == 'all') {
                            mform = paste0(contrast, ' ~ 1 + (1 | Subject) + (1 | fROI)')
                        } else {
                            mform = paste0(contrast, ' ~ 1')
                        }
                    } else { # mtype == 'abl
                        if (fROI == 'all') {
                            mform = paste0(contrast, ' ~ 0 + (1 | Subject) + (1 | fROI)')
                        } else {
                            mform = paste0(contrast, ' ~ 0')
                        }
                    }

                    if (fROI == 'all') {
                        m = lmer(mform, REML=F, data=df)
                    } else {
                        df_ = df[df$fROI == fROI,]
                        m = lm(mform, data=df_)
                    }

                    save(m, file=paste0(prefix_out, '.', contrast, '.', fROI, '.', mtype, '.lme.Rdata'))

                    sink(paste0(prefix_out, '.', contrast, '.', fROI, '.', mtype, '.lme.summary.txt'))
                    pr("Formula:", stdout())
                    pr(mform, stdout())
                    print(summary(m))
                    sink()
                }
            }
        }
    }
}

contrasts = c(
    'C6words_v_C1',
    'C6words_v_C2',
    'NLen6words_v_NLenC1',
    'NLen6words_v_NLenC2'
)


# Comparison between 6words and Nlength
for (parcel_set in c('evlab', 'PDD', 'RH', 'PDDanat')) {
    df_exp1 = read.table(paste0('contrasts/', parcel_set, '_nlength1_contrasts.csv'), header=TRUE, sep=',')
    df_exp2 = read.table(paste0('contrasts/', parcel_set, '_nlength2_contrasts.csv'), header=TRUE, sep=',')
    df_exp3 = read.table(paste0('contrasts/', parcel_set, '_6words_contrasts.csv'), header=TRUE, sep=',')
    prefix_out = paste0('glm/', parcel_set, '_betweengroups')

    for (contrast in contrasts) {
        if (startsWith(contrast, 'C')) {
            df_a = df_exp3[, c('Subject', 'fROI', 'C126_6words')]
        } else {
            df_a = df_exp3[, c('Subject', 'fROI', 'NLen6words')]
        }
        colnames(df_a)[3] = 'Effect'
        df_a$experiment = 1
        if (endsWith(contrast, '1')) {
            df_b = df_exp1
        } else {
            df_b = df_exp2
        }
        if (grepl('NLen', contrast)) {
            df_b = df_b[, c('Subject', 'fROI', 'NLenC')]
        } else {
            df_b = df_b[, c('Subject', 'fROI', 'C126')]
        }
        colnames(df_b)[3] = 'Effect'
        df_b$experiment = 0
        df = rbind(df_a, df_b)
        for (fROI in c('all', unique(df_exp3$fROI))) {
            for (mtype in c('full', 'abl')) {
                if (mtype == 'full') {
                    if (fROI == 'all') {
                        mform = 'Effect ~ experiment + (1 | Subject) + (1 | fROI)'
                    } else {
                        mform = 'Effect ~ experiment'
                    }
                } else { # mtype == 'abl
                    if (fROI == 'all') {
                        mform = 'Effect ~ 1 + (1 | Subject) + (1 | fROI)'
                    } else {
                        mform = 'Effect ~ 1'
                    }
                }

                if (fROI == 'all') {
                    m = lmer(mform, REML=F, data=df)
                } else {
                    df_ = df[df$fROI == fROI,]
                    m = lm(mform, data=df_)
                }

                save(m, file=paste0(prefix_out, '.', contrast, '.', fROI, '.', mtype, '.lme.Rdata'))

                sink(paste0(prefix_out, '.', contrast, '.', fROI, '.', mtype, '.lme.summary.txt'))
                pr("Formula:", stdout())
                pr(mform, stdout())
                print(summary(m))
                sink()
            }
        }
    }
}



# Pairwise fROI comparison

df = read.table(paste0('contrasts/evlab_nlength2_contrasts.csv'), header=TRUE, sep=',')
prefix_out = paste0('glm/evlab_nlength2')

# NLenC vs NLenJ in IFG vs PTL

contrast = 'NLenC_v_NLenJ_diff'
fROI_b = 'LPostTemp'
df_b = df[df$fROI == fROI_b, c('Subject', 'NLenC_v_NLenJ')]
df_b = df_b[order(df_b$Subject),]

for (fROI_a in c('LIFG', 'LIFGorb')) {
    df_a = df[df$fROI == fROI_a, c('Subject', 'NLenC_v_NLenJ')]
    df_a = df_a[order(df_a$Subject),]
    df_a[[contrast]] = df_a$NLenC_v_NLenJ - df_b$NLenC_v_NLenJ

    fROI = paste0(fROI_a, '_v_', fROI_b)
    for (mtype in c('full', 'abl')) {
        if (mtype == 'full') {
            mform = paste0(contrast, ' ~ 1')
        } else { # mtype == 'abl
            mform = paste0(contrast, ' ~ 0')
        }

        m = lm(mform, data=df_a)

        save(m, file=paste0(prefix_out, '.', contrast, '.', fROI, '.', mtype, '.lme.Rdata'))

        sink(paste0(prefix_out, '.', contrast, '.', fROI, '.', mtype, '.lme.summary.txt'))
        pr("Formula:", stdout())
        pr(mform, stdout())
        print(summary(m))
        sink()
    }
}


# NLenJ in AngG vs others

contrasts = c('NLenC', 'NLenJ', 'NLenC_v_NLenJ')
fROI_b = 'LAngG'

for (source_contrast in contrasts) {
    contrast = paste0(source_contrast, '_diff')
    for (fROI_a in c('LIFGorb', 'LIFG', 'LMFG', 'LAntTemp', 'LPostTemp')) {
        df_a = df[df$fROI == fROI_a, c('Subject', source_contrast)]
        df_a = df_a[order(df_a$Subject),]
        df_b = df[df$fROI == fROI_b, c('Subject', source_contrast)]
        df_b = df_b[order(df_b$Subject),]
        df_a[[contrast]] = df_a[[source_contrast]] - df_b[[source_contrast]]

        fROI = paste0(fROI_a, '_v_', fROI_b)
        for (mtype in c('full', 'abl')) {
            if (mtype == 'full') {
                mform = paste0(contrast, ' ~ 1')
            } else { # mtype == 'abl
                mform = paste0(contrast, ' ~ 0')
            }

            m = lm(mform, data=df_a)

            save(m, file=paste0(prefix_out, '.', contrast, '.', fROI, '.', mtype, '.lme.Rdata'))

            sink(paste0(prefix_out, '.', contrast, '.', fROI, '.', mtype, '.lme.summary.txt'))
            pr("Formula:", stdout())
            pr(mform, stdout())
            print(summary(m))
            sink()
        }
    }
}


# Laterality comparison

fROIs = c('IFGorb', 'IFG', 'MFG', 'AntTemp', 'PostTemp', 'AngG')

for (experiment in c(1, 2, 3)) {
    if (experiment < 3) {
        df_name = paste0('nlength', experiment)
    } else {
        df_name = '6words'
    }
    df1 = read.table(paste0('contrasts/evlab_', df_name, '_contrasts.csv'), header=TRUE, sep=',')
    df2 = read.table(paste0('contrasts/RH_', df_name, '_contrasts.csv'), header=TRUE, sep=',')
    df = rbind(df1, df2)

    prefix_out = paste0('glm/betweenhemispheres_', df_name)
    if (experiment == 1) {
        source_contrasts = c(
            'C',
            'NLenC',
            'S_v_W'
        )
    } else if (experiment == 2) {
        source_contrasts = c(
            'C',
            'J',
            'N',
            'C_v_J',
            'C_v_N',
            'NLenC',
            'NLenJ',
            'NLenN',
            'NLenC_v_NLenJ',
            'NLenC_v_NLenN',
            'NLenJ_v_NLenN',
            'S_v_W',
            'S_v_N',
            'J_v_W',
            'J_v_N',
            'S_v_W_v_J_v_N'
        )
    } else {
        source_contrasts = c(
            'C_WL24',
            'C_WL30',
            'C_6words',
            'NLen24WL',
            'NLen30WL',
            'NLen6words'
        )
    }

    for (source_contrast in source_contrasts) {
        contrast = paste0(source_contrast, '_diff')

        all = NULL

        for (fROI_src in fROIs) {
            fROI_a = paste0('L', fROI_src)
            fROI_b = paste0('R', fROI_src)
            df_a = df[df$fROI == fROI_a, c('Subject', source_contrast)]
            df_a = df_a[order(df_a$Subject),]
            df_b = df[df$fROI == fROI_b, c('Subject', source_contrast)]
            df_b = df_b[order(df_b$Subject),]
            df_a[[contrast]] = df_a[[source_contrast]] - df_b[[source_contrast]]
            df_a$fROI = paste0(fROI_src)
            if (is.null(all)) {
                all = df_a
            } else {
                all = rbind(all, df_a)
            }

            fROI = fROI_src
            for (mtype in c('full', 'abl')) {
                if (mtype == 'full') {
                    mform = paste0(contrast, ' ~ 1')
                } else { # mtype == 'abl
                    mform = paste0(contrast, ' ~ 0')
                }

                m = lm(mform, data=df_a)

                save(m, file=paste0(prefix_out, '.', contrast, '.', fROI, '.', mtype, '.lme.Rdata'))

                sink(paste0(prefix_out, '.', contrast, '.', fROI, '.', mtype, '.lme.summary.txt'))
                pr("Formula:", stdout())
                pr(mform, stdout())
                print(summary(m))
                sink()
            }
        }

        for (mtype in c('full', 'abl')) {
            if (mtype == 'full') {
                mform = paste0(contrast, ' ~ 1 + (1 | Subject) + (1 | fROI)')
            } else { # mtype == 'abl
                mform = paste0(contrast, ' ~ 0 + (1 | Subject) + (1 | fROI)')
            }

            m = lmer(mform, REML=F, data=all)

            save(m, file=paste0(prefix_out, '.', contrast, '.all.', mtype, '.lme.Rdata'))

            sink(paste0(prefix_out, '.', contrast, '.all.', mtype, '.lme.summary.txt'))
            pr("Formula:", stdout())
            pr(mform, stdout())
            print(summary(m))
            sink()
        }
    }
}

