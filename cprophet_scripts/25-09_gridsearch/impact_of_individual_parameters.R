require(gridExtra)
library(grid)
library(ggplot2)

source('read_data.R')
stats.file <- '/Volumes/rhafen/grid-search-runs/openms-complex-feat2/stats.tsv'
setwd('~/Dev/MAScripts/cprophet_scripts/25-09_gridsearch/')
stats <- readStatsData(stats.file)

stats.bak <- copy(stats)
target.vars <- c(
    'n_complexes_target_apmw_ok',
    'score_apmw',
    'protein_apex_spread_mean_target',
    'complexes_target_median_size'
)
param.col.idxs <- grep('^param__', colnames(stats.bak))
colnames(stats.bak) <- sub('param__', '', colnames(stats.bak))
param.names <- colnames(stats.bak)[param.col.idxs]

plots <- lapply(param.names, function(param.name) {
    target.plots <- lapply(target.vars, function(target.var) {
        formula <- as.formula(paste(target.var, '~', param.name))
        print(formula)
        x.mapping <- sprintf('as.factor(%s)', param.name)
        ggplot(stats.bak) + geom_boxplot(aes_string(x.mapping, target.var)) +
            xlab('')
    })
    do.call(arrangeGrob, c(target.plots,
                           ncol=length(target.plots),
                           top=param.name))
})
full.plot <- do.call(arrangeGrob, c(plots, nrow=length(plots)))
pdf('impact_individual_parameters_.pdf', 10, 80, onefile=T)
plot(full.plot)
dev.off()
