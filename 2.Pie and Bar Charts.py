#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 12:14:33 2021

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
import itertools 
from textwrap import wrap
import seaborn as sns


df = pd.read_csv("CellAge Senescence Genes.csv")
n_observations = df.shape[0]


# =====================================
# GRAPHS OF VARIABLES =================
# =====================================
#  METHOD
unique_methods = pd.unique(df["Method"]).tolist()

count_method = {}
for method in unique_methods:
    count_method[method] = df[df["Method"] == method].shape[0]

percentage_method = {}
for method in unique_methods:
    percentage_method[method] = (count_method[method] / n_observations) * 100

# # Check if the percentages add to 100
# sum_method_percentage = 0
# for method in percentage_method:
#     sum_method_percentage = sum_method_percentage + percentage_method[method]
# print(sum_method_percentage)

# Plot pie chart with the different methods
labels = unique_methods
sizes = percentage_method.values()
plt.pie(sizes, labels = labels, autopct='%1.1f%%')
plt.axis('equal')
plt.title("Pie chart of different methods")
plt.show()


# CELL TYPES
# Get unique values of column "Cell Types", sorted alphabetically
cell_types = df['Cell Types'].str.split(
    ',\s*', expand=True).stack().unique().tolist()
cell_types.sort()


# Count number of observations for each cell type
count_cell_types = {}
for category in cell_types:
    count_cell_types[category] = df["Cell Types"].str.count(category).sum()

# Sort the cell types by number of observations, in descending order
count_cell_types = {key: val for key, val in sorted(count_cell_types.items(), 
                                    key = lambda ele: ele[1], reverse = True)}

# Get first n items in dictionary 
n = 5
"""Chose 5 because after that the cell types accounts for less than 5% of total 
observations"""

top_cell_types = dict(itertools.islice(count_cell_types.items(), n))


# Top cell types in percentage
percentage_cell_types = {}
for cell_type in top_cell_types.keys():
    percentage_cell_types[cell_type] = (
        top_cell_types[cell_type] / n_observations) * 100

# Bar graph with the top cell types in the data
labels = list(top_cell_types.keys())
labels = [ '\n'.join(wrap(l, 10)) for l in labels ]
plt.title("Top 5 cell types in percentage")
plt.ylabel('Percentage of total observations')
# plt.grid(color='#95a5a6', linestyle='--', linewidth=2, axis='y', alpha=0.7)
plt.bar((labels), percentage_cell_types.values())
plt.show()


# CELL LINES 
# Get unique values of column "Cell Lines", sorted alphabetically
cell_lines_categories = df['Cell Lines'].str.split(
    ',\s*', expand=True).stack().unique().tolist()
cell_lines_categories.sort()

# Count number of observations for each cell line
count_cell_lines = {}
for category in cell_lines_categories:
    count_cell_lines[category] = df["Cell Lines"].str.count(category).sum()

# Sort the cell lines by number of observations, in descending order
count_cell_lines = {key: val for key, val in sorted(count_cell_lines.items(), 
                                    key = lambda ele: ele[1], reverse = True)}

# Get first n items in dictionary 
n = 6
"""Chose 6 because after that the cell types accounts for less than 5% of total 
observations"""

top_cell_lines = dict(itertools.islice(count_cell_lines.items(), n))


# Top cell lines in percentage terms
# n_total_unique_cell_lines = sum(count_cell_types.values())
percentage_cell_lines = {}
for cell_line in top_cell_lines.keys():
    percentage_cell_lines[cell_line] = (
        top_cell_lines[cell_line] / n_observations) * 100

# Bar graph with the top cell lines in the data
labels = list(top_cell_lines.keys())
labels = [ '\n'.join(wrap(l, 10)) for l in labels ]
plt.title("Top 6 cell lines in percentage")
plt.ylabel('Percentage of total observations')
# plt.grid(color='#95a5a6', linestyle='--', linewidth=2, axis='y', alpha=0.7)
plt.bar((labels), percentage_cell_lines.values())
plt.show()


# SENESCENSE TYPE
# Get unique values of column "Senescense Type", sorted alphabetically
senescence_types = df["Senescence Type"].unique().tolist()
senescence_types.sort()

# Count number of observations for each senescence type
count_senescence_types = {}
for category in senescence_types:
    count_senescence_types[category] = df["Senescence Type"].str.count(
        category).sum()

# Percentage number of observations per senescence type
percentage_senescense_types = {}
for senescense_type in count_senescence_types.keys():
    percentage_senescense_types[senescense_type] = (
        count_senescence_types[senescense_type] / n_observations) * 100
    
# Graph
x,y = 'Cancer Line?', 'Senescence Type'

df1 = pd.DataFrame(df.groupby(x)[y].value_counts())
df1.rename(columns={"Senescence Type":"percent"}, inplace = True)
df1["percent"] = df1["percent"] / n_observations
df1 = df1.mul(100)
df1 = df1.reset_index()

g = sns.catplot(x = x,y = 'percent', hue = y, kind = 'bar', data = df1, 
                aspect = 1.4, legend_out = False)

for p in g.ax.patches:
    txt = str(p.get_height().round(2)) + '%'
    txt_x = p.get_x() 
    txt_y = p.get_height()
    g.ax.text(txt_x,txt_y,txt)


# SENESCENSE EFFECT
# Get unique values of column "Senescense Effect", sorted alphabetically
senescence_effect = df["Senescence Effect"].unique().tolist()
senescence_effect.sort()

# Count number of observations for each senescence effect
count_senescence_effect = {}
for effect in senescence_effect:
    count_senescence_effect[effect] = df["Senescence Effect"].str.count(
        effect).sum()

# Percentage number of observations per senescence effect
percentage_senescense_effect = {}
for senescense_effect in count_senescence_effect.keys():
    percentage_senescense_effect[senescense_effect] = (
        count_senescence_effect[senescense_effect] / n_observations) * 100

# Graph
x,y = 'Cancer Line?', 'Senescence Effect'

df1 = pd.DataFrame(df.groupby(x)[y].value_counts())
df1.rename(columns={"Senescence Effect":"percent"}, inplace = True)
df1["percent"] = df1["percent"] / n_observations
df1 = df1.mul(100)
df1 = df1.reset_index()

h = sns.catplot(x=x,y='percent',hue=y,kind='bar',data=df1, aspect = 1, 
                legend_out = False)

for p in h.ax.patches:
    txt = str(p.get_height().round(2)) + '%'
    txt_x = p.get_x() 
    txt_y = p.get_height()
    h.ax.text(txt_x,txt_y,txt)

