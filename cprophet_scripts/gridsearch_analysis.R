# Remove zombie variables
# rm(list=ls())

# Load required packages and install in case they are not installed already.
required.packages <- c('ggplot2', 'rpart', 'rpart.plot', 'glmnet', 'reshape2',
                       'venneuler', 'data.table')
chooseCRANmirror(ind=82)
for (p in required.packages) {
    # Need to specify character.only=TRUE so that `require` will treat the
    # argument as a character variable.
    import.ok <- require(p, char=T)
    if (!import.ok) {
        install.packages(p)
        require(p, char=T)
    }
}

# Read output from gridsearch
stats.raw <- read.csv('/Volumes/rhafen/grid-search-runs/openms/stats.tsv', sep='\t')

## Compute rank based on number of identified target complexes
# Compute 'apmw ok' version of decoy rate
stats.raw$complex_decoy_rate_apmw_ok <- stats.raw$n_complexes_decoy_apmw_ok /
                                        stats.raw$n_complexes_target_apmw_ok


## Impute missing/zero/infinite values since they cause problems downstream
# Remove zero values from decoy rate
is.zero <- stats.raw$complex_decoy_rate_apmw_ok == 0
# Impute Inf/NaN (source: 0/0 == NaN, x/0 == Inf)
row.is.na <- is.na(stats.raw$complex_decoy_rate_apmw_ok)
is.inf <- is.infinite(stats.raw$complex_decoy_rate_apmw_ok)
stats.raw <- subset(stats.raw, !row.is.na & !is.inf & !is.zero)

stats.raw$diff_target_decoy_apwm_ok <- stats.raw$n_complexes_target_apmw_ok -
                                    stats.raw$n_complexes_decoy_apmw_ok
stats.raw$score_apmw <- 1 / stats.raw$complex_decoy_rate_apmw_ok * stats.raw$diff_target_decoy_apwm_ok
stats.raw$score <- 1 / stats.raw$complex_decoy_rate * stats.raw$diff_target_decoy

par.backup <- par(mfrow=c(3,1))

## Set the target variables
target.var <- 'score_apmw'
# target.var <- 'score'

## Order the rows of the data set according to the target variable
target.var.order <- order(stats.raw[[target.var]], decreasing=TRUE)
stats.raw <- stats.raw[target.var.order, ]

stats.openms <- subset(stats.raw, param__protfeatfinder == 'openms')
stats.spline <- subset(stats.raw, param__protfeatfinder == 'spline')

## Set the stats dataset that is used in downstream analysis
stats <- stats.openms

## Plot
if (target.var == 'score') {
    plot(stats$n_complexes_decoy, col='blue')
    points(stats[[target.var]])
    points(stats$complex_decoy_rate, col='green', pch=3)
} else {
    plot(stats[[target.var]])
    points(stats$n_complexes_decoy_apmw_ok, col='blue')
    points(stats$complex_decoy_rate_apmw_ok, col='green', pch=3)
}

## Variable selection
# Create training set
getParams <- function(df) {
  is.param.col <- grepl('^param__', colnames(df))
  subset(df, select=is.param.col)
}
is.target.col <- colnames(stats) == target.var
is.param.col <- grepl('^param__', colnames(stats))
trainset <- stats[, is.target.col | is.param.col]
colnames(trainset) <- sub('param__', '', colnames(trainset))

# Decision tree
fit <- rpart(as.formula(paste(target.var, '~ .')), data=trainset)
rpart.plot(fit, varlen=0)

# LASSO
trainset.lasso <- trainset
X <- trainset.lasso[, colnames(trainset) != target.var]
y <- trainset.lasso[[target.var]]
is.fac <- sapply(X, is.factor)
X[is.fac] <- data.frame(do.call(cbind, lapply(X[, is.fac], as.integer)))
non.changing <- sapply(X, function(col) sd(col) == 0)
X <- as.matrix(X[, !non.changing])
X <- scale(X)
lasso.cv <- cv.glmnet(X, y)
coef(lasso.cv)

# Read corum identifiers
corum <- read.csv('~/Dev/MAScripts/data/corum_complex_protein_assoc.tsv', sep='\t',
                  stringsAsFactors=FALSE)
corum.complexes <- unique(corum$complex_id)

plotVennDiagram <- function(stats, n.topx=3, with.apmw=FALSE) {
    if (with.apmw) {
        ident.complexes <- strsplit(as.character(stats$target_complex_ids_apwm_ok), ',')
        targetvar.order <- order(stats$score_apmw, decreasing=T)
    } else {
        ident.complexes <- strsplit(as.character(stats$target_complex_ids), ',')
        targetvar.order <- order(stats$score, decreasing=T)
    }
    ## Plot overlap between detected corum complexes
    # Select the detected complex sets of the X parameter sets that have the highest
    # difference between target and decoys.
    stats.topx <- stats[targetvar.order[1:n.topx], ]
    ident.complexes.topx <- ident.complexes[targetvar.order[1:n.topx]]
    # Add the iteration number as the list's names (will later be used as the set name)
    names(ident.complexes.topx) <- as.character(targetvar.order[1:n.topx]-1)
    # Add all the corum complexes as a last element
    ident.complexes.topx$CORUM <- corum.complexes
    # Convert the list to a long format df
    dat <- melt(ident.complexes.topx, value.name='elements', variable.name='sets')
    v <- venneuler(dat)
    plot(v)
}


## Plot overlap between detected corum complexes
# Select the detected complex sets of the X parameter sets that have the highest
# difference between target and decoys.
plotVennDiagram(stats, n.topx=3, with.apmw=target.var == 'score_apmw')

par(par.backup)
