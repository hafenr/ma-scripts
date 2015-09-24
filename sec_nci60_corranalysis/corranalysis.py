import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

from rpy2 import robjects as ro
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.lib import grid
from rpy2.robjects.lib import ggplot2 as gg
from rpy2.robjects.packages import importr

# Execute the scripts that read the nci60 correlations and compute
# the average protein trace correlation for each protein complex.
from sec_data import complex_corr_sec
from nci60_data import complex_corr_nci60

# Check what complexes were identified in both experiments
nci60_complexes = set(complex_corr_nci60.index)
sec_complexes = set(complex_corr_sec.index)
complexes_in_both = nci60_complexes.intersection(sec_complexes)
print 'NUM COMPLEXES IN BOTH:\t%d' % len(complexes_in_both)

# Subset the data so that we only look at complexes that are in both experiments
corr_sec = complex_corr_sec[complexes_in_both]
corr_nci60 = complex_corr_nci60[complexes_in_both]

# How do they correlate?
total_correl = np.corrcoef(corr_sec, corr_nci60)[0, 1]
print 'TOTAL CORRELATION:\t%f' % total_correl

# Necessary so that rpy2 can automatically convert pandas dataframes into
# R dataframes.
pandas2ri.activate()

# Produce a linear fit
x = corr_nci60
y = corr_sec
xmin = np.min(x)
xmax = np.max(x)
xs = np.linspace(xmin, xmax, num=100).reshape(100, 1)
lm = LinearRegression()
# The training data for scikit models must be in matrix
# form, i.e. columns == features, rows == observations.
# For this we need to reshape the 1-dimensional arrays.
X = corr_nci60.reshape(len(x), 1)
y = corr_sec
lm.fit(X, y)
y_pred = lm.predict(xs)

# Plot the data using the R-bridge rpy and ggplot
p = gg.ggplot(pd.DataFrame())
p += gg.geom_point(
    gg.aes_string(x='r_nci60', y='r_sec'),
    data=pd.DataFrame({
        'r_nci60': corr_nci60,
        'r_sec': corr_sec
    })
)
p += gg.geom_line(
    gg.aes_string(x='x', y='y'),
    data=pd.DataFrame({
        'x': xs.reshape(-1),
        'y': y_pred
    }),
    color='red'
)
p.plot()
