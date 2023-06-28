#!/usr/bin/env Rscript

pr = function(x, buffer=NULL) {
    if (is.null(buffer)) {
        buffer = stderr()
    }
    cat(paste0(x, '\n'), file=buffer, append=TRUE)
}

library(lme4)


for (parcel_set in c('betweenhemispheres', 'evlab', 'PDD', 'RH')) {
    for (experiment in c(1, 2, 3, 4)) {
        if (experiment < 3) {
            df_name = paste0('nlength', experiment)
        } else if (experiment == 3) {
            df_name = '6words'
        } else {
            df_name = 'betweengroups'
        }
        prefix_in = paste0('glm/', parcel_set, '_', df_name)
        prefix_out = paste0('lrt/', parcel_set, '_', df_name)

        if (!dir.exists('lrt')) {
            dir.create('lrt', recursive = TRUE)
        }

        models = list.files(path='glm', pattern=paste0(parcel_set, '_', df_name, '.*.Rdata'), full.names=TRUE, recursive=FALSE)

        for (mpath in models) {
            if (grepl('full', mpath)) {
                path_parts = strsplit(mpath, '.', fixed=TRUE)[[1]]
                fROI = path_parts[[length(path_parts) - 3]]
                contrast = path_parts[[length(path_parts) - 4]]

                m_full_path = mpath
                m_abl_path = gsub('full', 'abl', mpath)

                m_full = get(load(m_full_path))
                m_abl = get(load(m_abl_path))
                lrt = anova(m_full, m_abl)

                sum_path = paste0(prefix_out, '.', contrast, '.', fROI, '.lrt.summary.txt')

                sink(sum_path)
                pr('==================\nLikelihood ratio test\n', stdout())
                pr(paste0('Experiment: ', df_name), stdout())
                pr(paste0('Variable:   \n', contrast), stdout())
                pr(paste0('fROI:       \n', fROI), stdout())
                print(lrt)
                pr('------------------\nFull model\n', stdout())
                pr(paste0('Path:       ', m_full_path), stdout())
                print(summary(m_full))
                pr('\n\n', stdout())
                pr('------------------\nAblated model\n', stdout())
                pr(paste0('Path:       ', m_abl_path), stdout())
                print(summary(m_abl))
                sink()
            }
        }
    }
}

