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


import plotly.express as px
import pandas as pd


desired_width = 320
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', desired_width)


df = pd.read_csv("CellAge Senescence Genes.csv")


# PARALLEL CATEGORY PLOT
l = []
for rowIndex, row in df.iterrows():
    if row["Senescence Effect"] == "Induces":
        l.append(1)
    elif row["Senescence Effect"] == "Inhibits":
        l.append(-1)
    else:
        l.append(0)
   
df1 = df.copy()
df1["Senescense Effect Number"] = l

# Graph
categories = ['Method', 'Cancer Line?', 'Senescence Type']
fig = px.parallel_categories(df1, dimensions = categories, 
                             color = "Senescense Effect Number")
fig.show()
print("Explanation of the graph: each vertical column represents a", 
      "category or column of the spreadsheet. In this case I chose the", 
      "'Method', 'Cancer Line?' and 'Senescence Type' columns/variables. The",
      "column to the right of the graph attributes a color to a number,", 
      "which then corresponds to a senescense effect. So, 'Induces' = 1,",
      "'Inhibits = -1' and 'Unclear' = 0. Each band in the graph shows a", 
      "permutation of the different categories/variables.")
print("\nUnfortunately the categories of each column are not entirely",
      "visible but there are no indications on how to change this in the",
      "documentation. The descriptions correspoding to the 'Senescence Effect",
      "Number' are also not visible for the same reason. But they are,", 
      "again: Induces = 1, Inhibits = -1, Unclear = 0.")
print("\nAs an example,", 
      "the graph shows that there are 59 cases (the biggest band), whose", 
      "method was overexpression, without cancer line, replicative with", 
      "regards to senescence type and in which there was an inducing",
      "senescense effect.")


# TREEMAP
categories = ['Method', 'Cancer Line?', 'Senescence Type', "Senescence Effect"]

df2 = df.groupby(categories).sum()
df2 = df2.assign(count = df.groupby(categories)["Method"].count())
df2.reset_index(drop = False, inplace = True)

l = []
for rowIndex, row in df2.iterrows():
    if row["Senescence Effect"] == "Induces":
        l.append(1)
    elif row["Senescence Effect"] == "Inhibits":
        l.append(-1)
    else:
        l.append(0)
   
df3 = df2.copy()
df3["Senescense Effect Number"] = l

# Graph
fig = px.treemap(df3, path = categories, values = "count", 
                  color = 'Senescence Effect')
fig.show()
print("\n\nAttempt at producing a Treemap. It's probably not very",
      "adequate for this particular case because many combinations of",
      "variables only have a few observation.")

