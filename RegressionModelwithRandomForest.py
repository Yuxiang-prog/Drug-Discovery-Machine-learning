import pandas as pd
import seaborn as sns
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
import matplotlib.pyplot as plt

df = pd.read_csv('acetylcholinesterase_06_bioactivity_data_3class_pIC50_pubchem_fp.csv')
X = df.drop('pIC50', axis=1)
# creates x-variable matrix -- only contains pubchem fingerprints
Y = df.pIC50
# creates y-variable matrix

from sklearn.feature_selection import VarianceThreshold

selection = VarianceThreshold(threshold=(.8 * (1 - .8)))
X = selection.fit_transform(X)
# remove low variance features

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2)
# split the data ina 80-20 fashion

# Building a regression model using random forest
model = RandomForestRegressor(n_estimators=100)
np.random.seed(100)
model.fit(X_train, Y_train)
r2 = model.score(X_test, Y_test)
print(r2)

Y_pred = model.predict(X_test)

sns.set(color_codes=True)
sns.set_style("white")

ax = sns.regplot(Y_test, Y_pred, scatter_kws={'alpha': 0.4})
ax.set_xlabel('Experimental pIC50', fontsize='large', fontweight='bold')
ax.set_ylabel('Predicted pIC50', fontsize='large', fontweight='bold')
ax.set_xlim(0, 12)
ax.set_ylim(0, 12)
ax.figure.set_size_inches(5, 5)
plt.show()
# this segment of code creates a scatter plot of the experimental vs predicted
# pIC50 values
# A higher r^2 value indicates that the model has found the correct characteristics
# for a more potent drug
