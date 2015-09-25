required.packages <- c('ggplot2', 'rpart', 'rpart.plot', 'glmnet', 'reshape2',
                       'venneuler', 'data.table')
sapply(required.packages, function(p) require(p, char=T))


setwd('~/Dev/MAScripts/cprophet_scripts/25-09_gridsearch/')


# 1. Plot the ranked parameter sets plot
source('read_data.R')
source('plot_ranked_param_sets.R')

stats <- readStatsData('/Volumes/rhafen/grid-search-runs/openms/stats.tsv')
plotRankedParameterSets(stats, 'score_apmw')

# Decision tree
source('variable_selection.R')
fitDecisionTree(stats, 'score_apmw')

