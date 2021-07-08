#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 10:42:20 2021

@author: Diogo
"""

"""
CLEAR THE CONSOLE AND REMOVE ALL VARIABLES PRESENT ON THE NAMESPACE. THIS IS 
USEFUL TO PREVENT PYTHON FROM CONSUMING MORE RAM EACH TIME I RUN THE CODE.
"""
try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass


import pandas as pd
from sklearn.preprocessing import OneHotEncoder, LabelBinarizer
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, cross_val_predict, cross_validate
from sklearn.metrics import accuracy_score, precision_score, roc_curve, auc, make_scorer, roc_auc_score
from sklearn.linear_model import LogisticRegression
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv("CellAge Senescence Genes.csv")

# DROP OBSERVATIONS WHERE "SENESCENCE EFFECT" = "UNCLEAR"
for rowIndex, row in df.iterrows():
    if row["Senescence Effect"] == "Unclear":
        df.drop(rowIndex, axis = 0, inplace = True)

# =====================================
# FIT LOGISTIC REGRESSION =============
# =====================================
# ENCONDE CATEGORICAL INDEPENDENT VARIABLES (CONVERT THEM INTO NUMERIC VARIABLES)
ohe = OneHotEncoder()
df2 = df.copy()
df2.drop(["Entrez Id", "Gene Symbol", "Senescence Effect"], axis = 1,
         inplace = True)
ohe.fit(df2)
X = ohe.transform(df2).toarray()

"""ENCODE CATEGORICAL DEPENDENT VARIABLE (SENESCENCE EFFECT) (CONVERT IT TO A
NUMERIC VARIABLES)"""
lb = LabelBinarizer()
lb.fit(df["Senescence Effect"])
y = np.ravel(lb.transform(df["Senescence Effect"])) 
# ravel() returns a 1D array of the input

# SPLIT DATA INTO TRAIN AND TEST SAMPLES
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1)


"""FITTING LOGISTIC REGRESSION TO ENCODED DATA AND CALCULATING THE MODEL'S 
ACCURACY"""
lr = LogisticRegression().fit(X_train, y_train)

# How to decode the data
# ohe.inverse_transform(X[0, :])
# lb.inverse_transform(y[0])

# # Check if the model is returning probabilities
# lr.predict_proba(X[0, :].reshape(1,-1)) 
# # Returns the probabilities for the 2 classes, adding to 1.

# PREDICTIONS AND ACCURACY
"""Accuracy here is defined by the ratio of correct predictions to total 
predictions. Technically it is the Jaccard similarity coefficient score."""  
predictions = lr.predict(X_test)
print("Relative (%) accuracy of the model when using one test sample:", 
      accuracy_score(predictions, y_test) * 100)

scores = cross_val_score(lr, X, y, cv = 5)
print("Relative (%) average accuracy across 5 different samples:", 
      scores.mean() * 100)
print("% Minimum and Maximum accuracy scores:", scores.min() * 100, 
      scores.max() * 100)
print("% Standard deviation:", scores.std() * 100)

"""The precision is the ratio tp / (tp + fp) where tp is the number of true 
positives and fp the number of false positives. Intuitively, precision is the
ability of the classifier not to label as positive a sample that is negative. 
The best value is 100% and the worst value is 0."""
print("\n% Precision Score:", precision_score(y_true = y_test, 
                                               y_pred = predictions) * 100)
scores2 = cross_val_score(lr, X, y, cv = 5, scoring = "precision")
print("% Cross Validated Precision Score:", scores2.mean() * 100)

# CONFUSION MATRIX
"""Each row of the matrix represents the instances in an actual class while
each column represents the instances in a predicted class. In our case 0 
represents an inducing senescence effect and 1 represents an inhibiting effect"""
confusion_matrix = pd.crosstab(y_test, predictions, rownames = ['Actual'], 
                                colnames = ['Predicted'])
sns.heatmap(confusion_matrix, annot=True)
plt.show()


# =============================================================================
# FIT ANOTHER LOGISTIC REGRESSION WITH LESS INDEPENDENT VARIABLES =============
# =============================================================================
# ENCONDE CATEGORICAL INDEPENDENT VARIABLES (CONVERT THEM INTO NUMERIC VARIABLES)
ohe = OneHotEncoder()
df3 = df.copy()
df3.drop(["Entrez Id", "Gene Symbol", "Cell Types", "Cell Lines", 
          "Cancer Line?", "Senescence Effect"], axis = 1, inplace = True)
ohe.fit(df3)
X2 = ohe.transform(df3).toarray()

# SPLIT DATA INTO TRAIN AND TEST SAMPLES
X_train, X_test, y_train, y_test = train_test_split(X2, y, test_size = 0.1)

# FITTING
lr = LogisticRegression().fit(X_train, y_train)


# PREDICTIONS, ACCURACY AND PRECISION
"""Accuracy here is defined by the ratio of correct predictions to total 
predictions. Technically it is the Jaccard similarity coefficient score.""" 
predictions = lr.predict(X_test)
print("\nRelative (%) accuracy of the model when using one test sample:", 
      accuracy_score(predictions, y_test) * 100)

scores = cross_val_score(lr, X2, y, cv = 5)
print("Relative (%) average accuracy across 5 different samples:", 
      scores.mean() * 100)
print("% Minimum and Maximum accuracy scores:", scores.min() * 100, 
      scores.max() * 100)
print("% Standard deviation:", scores.std() * 100)

print("\n% Precision Score:", precision_score(y_true = y_test, 
                                               y_pred = predictions) * 100)
scores2 = cross_val_score(lr, X2, y, cv = 5, scoring = "precision")
print("% Cross Validated Precision Score:", scores2.mean() * 100)

print("\nConclusion: removing some of the independent variables slightly", 
      "improved the model's accuracy and precision")

# CONFUSION MATRIX
"""Each row of the matrix represents the instances in an actual class while
each column represents the instances in a predicted class. In our case 0 
represents an inducing senescence effect and 1 represents an inhibiting 
effect."""
confusion_matrix = pd.crosstab(y_test, predictions, rownames = ['Actual'], 
                                colnames = ['Predicted'])
sns.heatmap(confusion_matrix, annot=True)
plt.show()


# COMPUTE RECEIVER OPERATING CHARACTERISTIC (ROC) and Area Under the Curve (AUC)
"""ROC curves typically feature true positive rate on the Y axis, and false 
positive rate on the X axis. This means that the top left corner of the plot 
is the “ideal” point - a false positive rate of zero, and a true positive rate
 of one. This is not very realistic, but it does mean that a larger area under
 the curve (AUC) is usually better."""
"""The “steepness” of ROC curves is also important, since it is ideal to 
maximize the true positive rate while minimizing the false positive rate."""
"""ROC curves are typically used in binary classification to study the output 
of a classifier."""
positive_class = 1 # Senescence effect = Inhibits = 1

y_score2 = cross_val_predict(lr, X, y, method='predict_proba')

fpr = dict()
tpr = dict()
roc_auc = dict()
fpr[positive_class], tpr[positive_class], _ = roc_curve(y, y_score2[:,
                                                        positive_class])
roc_auc[positive_class] = auc(fpr[positive_class], tpr[positive_class])
print("\nAUC:", roc_auc[positive_class])

myscore = make_scorer(roc_auc_score, needs_proba=True)

my_value = cross_validate(lr, X, y, scoring = myscore)
print("Cross-validated AUC:", np.mean(my_value['test_score'].tolist()))


# PLOT ROC CURVE FOR THE "INHIBITS" CLASSIFICATION
plt.figure()
lw = 2
plt.plot(fpr[1], tpr[1], color = 'darkorange',
         lw = lw, label = 'ROC curve (area = %0.2f)' % roc_auc[1])
plt.plot([0, 1], [0, 1], color = 'navy', lw = lw, linestyle = '--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.legend(loc="lower right")
plt.show()

