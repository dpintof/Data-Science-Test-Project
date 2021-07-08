#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 13:18:30 2021

@author: Diogo
"""

import json
import copy

def read_ipynb(notebook_path):
    with open(notebook_path, 'r', encoding='utf-8') as f:
        return json.load(f)
    
one = read_ipynb('1.Data Exploration.ipynb')
two = read_ipynb('2.Pie and Bar Charts.ipynb')
three = read_ipynb('3.Parallel Category Plot and Treemap.ipynb')
four = read_ipynb('4.Logistic Regression.ipynb')
five = read_ipynb('5.SVM Model.ipynb')
six = read_ipynb('6.Random Forest.ipynb')

final_notebook = copy.deepcopy(one)
final_notebook['cells'] = one['cells'] + two['cells'] + three['cells'] + four[
    'cells'] + five['cells'] + six['cells']

def write_ipynb(notebook, notebook_path):
    with open(notebook_path, 'w', encoding='utf-8') as f:
        json.dump(notebook, f)

# Saving the resulting notebook
write_ipynb(final_notebook, 'Ageless Jupyter Notebook.ipynb')

