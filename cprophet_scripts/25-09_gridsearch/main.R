required.packages <- c('ggplot2', 'rpart', 'rpart.plot', 'glmnet', 'reshape2',
                       'venneuler', 'data.table')
sapply(required.packages, function(p) require(p, char=T))

getParams <- function(df) {
  df <- subset(df, select=grepl('^param__', colnames(df)))
  colnames(df) <- sub('param__', '', colnames(df))
  df
}

source('read_data.R')
setwd('~/Dev/MAScripts/cprophet_scripts/25-09_gridsearch/')
stats.file <- '/Volumes/rhafen/grid-search-runs/openms-complex-feat2/stats.tsv'
stats <- readStatsData(stats.file)

# 1. Plot the ranked parameter sets plot
source('plot_ranked_param_sets.R')
plotRankedParameterSets(stats, 'score_apmw', xlim=c(0, 200))

# 2. Plot overlap between detected corum complexes
# Select the detected complex sets of the X parameter sets that have the highest
# difference between target and decoys.
source('complex_id_intersections.R')

# 3. How do the top comparing parameter sets compare in terms of their parameter values?
source('differing_parameters.R')

# 4. Parameter sets of top scoring runs
stats.dt <- as.data.table(stats)
params <- getParams(stats.dt)
n.topx <- 5
params.top <- params[1:n.topx, ]
print(params.top)
