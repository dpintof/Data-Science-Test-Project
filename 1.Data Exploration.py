#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 08:19:50 2021

@author: Diogo
"""

"""
Clear the console and remove all variables present on the namespace. This is 
useful to prevent Python from consuming more RAM each time I run the code.
"""
try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency
import numpy as np


df = pd.read_csv("CellAge Senescence Genes.csv")

# =====================================
# DATA EXPLORATION ====================
# =====================================

# CHECK FOR UNIQUE VALUES IN THE FEATURES/COLUMNS
# If a column has unique values it means there is not point in making a graph
print("Are there null values in the data?", df.isnull().values.any())
print("Are the strings in the Gene Symbol Column unique? ", 
      df["Gene Symbol"].is_unique)
print("Are the strings in the Cell Types Column unique? ", 
      df["Cell Types"].is_unique)
print("Are the strings in the Cell Lines Column unique? ", 
      df["Cell Lines"].is_unique)


# CORRELATIONS
"""Some variables seem to be correlated. Let's test that assumption using a 
Chi-square test"""
# Assumption(H0): The two columns are not correlated
cross_tab_result = pd.crosstab(index = df["Method"], columns = df[
    "Cell Types"])
chi_sq_result = chi2_contingency(cross_tab_result)

"""The P-Value for the chi-square test for independence is the probability of 
reaching a value (the chi-square statistic) different than expected, assuming 
that the two variables are not related (which is claimed by the null 
hypothesis). The smaller the p-value, the less likely it is that one would
observe a value so different than expected, assuming H0 is true. In other 
words, the smaller the p-value, the less likely it is that we are rejecting H0 
by pure chance."""

# If P-Value>0.05 then we fail to reject the null hypothesis (H0).
print('\nThe P-Value of the Chi-square test is:', chi_sq_result[1])
print("P-value > 0,05 so we fail to reject the null hypothesis that the",
      "variables 'Method' and 'Cell Types' are not correlated.")


# CORRELATIONS BEYOND 2 VARIABLES
df3 = df.copy()
df3.drop(["Entrez Id", "Gene Symbol", "Senescence Effect"], axis = 1, 
         inplace = True)

factors_paired = [(i,j) for i in df3.columns.values for j in df3.columns.values]
chi_sq_results, p_values =[], []

for f in factors_paired:
    if f[0] != f[1]:
        chitest = chi2_contingency(pd.crosstab(df3[f[0]], df3[f[1]]))   
        chi_sq_results.append(chitest[0])
        p_values.append(chitest[1])
    else:      # for same factor pair
        chi_sq_results.append(0)
        p_values.append(0)

p_values = np.array(p_values).reshape((5,5)) # shape it as a matrix
p_values = pd.DataFrame(p_values, index = df3.columns.values, 
                        columns=df3.columns.values)


# HEATMAP OF CORRELATIONS
sns.heatmap(p_values, annot=True)
print("\nFrom the heatmap of p-values we can see that only 'Cell Types' and", 
      "'Senescence Type', and 'Method' and 'Cell Types' are not correlated.", 
      "This implies that maybe we could drop, for example, the variables", 
      "'Cell Types', 'Cell Lines' and 'Cancer Lines' when creating models.") 
plt.show()

