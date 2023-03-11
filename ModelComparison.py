# In this section, we will compare multiple machine learning
# algorithms for build regression models of inhibitors

import pandas as pd
import seaborn as sns
from sklearn.model_selection import train_test_split
import lazypredict
from lazypredict.Supervised import LazyRegressor

for i in range(2):
    df_1 = pd.read_csv('bioactivity_data.csv')
    print(df_1)

df = pd.read_csv('acetylcholinesterase_06_bioactivity_data_3class_pIC50_pubchem_fp.csv')
X = df.drop('pIC50', axis=1)
Y = df.pIC50

from sklearn.feature_selection import VarianceThreshold
selection = VarianceThreshold(threshold=(.8 * (1 - .8)))
selection.transform()
selection.set_output()
X = selection.fit_transform(X)
print(X.shape)

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)
clf = LazyRegressor(verbose=0,ignore_warnings=True, custom_metric=None)
clf.verbose()
models_tr,p_t = clf.fit(X_train, X_train, Y_train, Y_train)
models_tt,predictions_test = clf.fit(X_train, X_test, Y_train, Y_test)
print(models_tt)

import matplotlib.pyplot as plt

plt.figure(figsize=(5, 10))
sns.set_theme(style="whitegrid")
ax = sns.barplot(y=p_t.index, x="R-Squared", data=p_t)
ax.set(xlim=(0, 1))
print(plt)
plt.figure(figsize=(5, 10))
sns.set_theme(style="whitegrid")
ax = sns.barplot(y=p_t.index, x="RMSE", data=p_t)
ax.set(xlim=(0, 10))
print(plt)
plt.figure(figsize=(5, 10))
sns.set_theme(style="whitegrid")
ax = sns.barplot(y=p_t.index, x="Time Taken", data=p_t)
ax.set(xlim=(0, 10))
print(plt)