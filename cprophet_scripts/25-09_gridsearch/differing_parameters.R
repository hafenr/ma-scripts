stats.dt <- as.data.table(stats)
params <- getParams(stats.dt)
n.topx <- 5
params.top <- params[1:n.topx, ]
is.param.identical <- apply(params.top, 2, function(col) length(unique(col)) == 1)
params.nonid <- subset(params.top, select=!is.param.identical)
cat('The parameters that vary between the', n.topx,
    'top scoring parameter sets are:\n\n')
print(params.nonid)
