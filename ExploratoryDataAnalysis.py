import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski


def lipinski(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:

        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])

        if (i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i + 1

    columnNames = ["MW", "LogP", "NumHDonors", "NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors
def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i * (10 ** -9)  # Converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)

    return x
def norm_value(input):
    norm = []

    for i in input['standard_value']:
        if i > 100000000:
            i = 100000000
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input.drop('standard_value', 1)

    return x

df = pd.read_csv('processeddata.csv')
df_no_smiles = df.drop(columns='canonical_smiles')
# dropping nonessential data

smiles = []
for i in df.canonical_smiles.tolist():
  cpd = str(i).split('.')
  longest = max(cpd, key=len)
  print(longest)
  smiles.append(longest)

smiles = pd.Series(smiles, name='canonical_smiles')
df_clean_smiles = pd.concat([df_no_smiles,smiles], axis=1)
df_clean_smiles

# Lipinski's Rule
# Evaluates the drug-likeness of a specific molecule

df_lipinski = lipinski(df_clean_smiles.canonical_smiles)
print(df_lipinski)

df_combined = pd.concat([df,df_lipinski], axis=1)
# combine the 2 DataFrames

df_combined.standard_value.describe()

df_norm = norm_value(df_combined)
df_final = pIC50(df_norm)
print(smiles)
print(df_final.pIC50.describe())

df_final.to_csv('acetylcholinesterase_04_bioactivity_data_3class_pIC50.csv')
df_2class = df_final[df_final['class'] != 'intermediate']
# removing the intermediate class from the dataset

df_2class.to_csv('acetylcholinesterase_05_bioactivity_data_2class_pIC50.csv')

# Performing chemical space analysis
import seaborn as sns
sns.set(style='ticks')
import matplotlib.pyplot as mat

mat.figure(figsize=(5.5, 5.5))

sns.countplot(x='class', data=df_2class, edgecolor='black')

mat.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
mat.ylabel('Frequency', fontsize=14, fontweight='bold')
mat.savefig('plot_bioactivity_class.pdf')
# just offers a visual representation of which inhibitors are labeled as "active" and "inactive"

mat.figure(figsize=(5.5, 5.5))
sns.scatterplot(x='MW', y='LogP', data=df_2class, hue='class', size='pIC50', edgecolor='black', alpha=0.7)

mat.xlabel('MW', fontsize=14, fontweight='bold')
mat.ylabel('LogP', fontsize=14, fontweight='bold')
mat.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
mat.savefig('plot_MW_vs_LogP.pdf')
# offers a scatterplot representation to begin training the model
# is there a relationship between pIC50 value, logP and if the drug inhibitor is active/inactive?

mat.figure(figsize=(5.5, 5.5))
sns.boxplot(x = 'class', y = 'pIC50', data = df_2class)
mat.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
mat.ylabel('pIC50 value', fontsize=14, fontweight='bold')
mat.savefig('plot_ic50.pdf')
print(smiles)

def mannwhitney(descriptor, verbose=False):
    # https://machinelearningmastery.com/nonparametric-statistical-significance-tests-in-python/
    from numpy.random import seed
    from numpy.random import randn
    from scipy.stats import mannwhitneyu

    # seed the random number generator
    seed(1)

    # actives and inactives
    selection = [descriptor, 'class']
    df = df_2class[selection]
    active = df[df['class'] == 'active']
    active = active[descriptor]

    selection = [descriptor, 'class']
    df = df_2class[selection]
    inactive = df[df['class'] == 'inactive']
    inactive = inactive[descriptor]

    # compare samples
    stat, p = mannwhitneyu(active, inactive)
    # print('Statistics=%.3f, p=%.3f' % (stat, p))

    # interpret
    alpha = 0.05
    if p > alpha:
        interpretation = 'Same distribution (fail to reject H0)'
    else:
        interpretation = 'Different distribution (reject H0)'

    results = pd.DataFrame({'Descriptor': descriptor,
                            'Statistics': stat,
                            'p': p,
                            'alpha': alpha,
                            'Interpretation': interpretation}, index=[0])
    filename = 'mannwhitneyu_' + descriptor + '.csv'
    results.to_csv(filename)

    return results

# the code segment above uses a machine learing model to test if there is a statistical difference between active and
# inactive inhibitors in term of the pIC50 variable -- remember the IC50/PIC50 is nothing more than a value that
# indicates the potency of a drug -- "how much drug is needed to inhibit a biological process by half"


mannwhitney('pIC50')
# since the p-value is extremely low, we can reject the null hypothesis

# do the same process for all the other gathered values
mat.figure(figsize=(5.5, 5.5))
sns.boxplot(x = 'class', y = 'MW', data = df_2class)
print(smiles)
mat.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
mat.ylabel('MW', fontsize=14, fontweight='bold')
mat.savefig('plot_MW.pdf')
mannwhitney('MW')

mat.figure(figsize=(5.5, 5.5))
sns.boxplot(x = 'class', y = 'LogP', data = df_2class)
print(smiles)
mat.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
mat.ylabel('LogP', fontsize=14, fontweight='bold')
mat.savefig('plot_LogP.pdf')
mannwhitney('LogP')

mat.figure(figsize=(5.5, 5.5))
sns.boxplot(x = 'class', y = 'NumHDonors', data = df_2class)
mat.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
mat.ylabel('NumHDonors', fontsize=14, fontweight='bold')
mat.savefig('plot_NumHDonors.pdf')
print(smiles)
mannwhitney('NumHDonors')


mat.figure(figsize=(5.5, 5.5))
sns.boxplot(x = 'class', y = 'NumHAcceptors', data = df_2class)
mat.xlabel('Bioactivity class', fontsize=14, fontweight='bold')
mat.ylabel('NumHAcceptors', fontsize=14, fontweight='bold')
print(smiles)
mat.savefig('plot_NumHAcceptors.pdf')
print(mannwhitney('NumHAcceptors'))

# all of the above distributions displayed "significant statistical difference"
