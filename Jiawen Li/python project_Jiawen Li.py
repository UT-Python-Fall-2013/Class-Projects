# -*- coding: utf-8 -*-
"""
Least square fit of an arbitrary nonlinear regression model to a given dataset in CSV format.

Created on Wed Dec 04 08:44:35 2013

@author: Jiawen Li
"""
import csv
import numpy as np
import pylab
from lmfit import minimize, Parameters, Parameter, report_fit

# some settings for the figure
font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 10}

pylab.rc('font', **font)

# Step 1. Get input file and parse the data
filename = "testdata.csv"

colors = "rgbkcmy" #different colors for plotting purposes

# =====================READ THE CSV INTO columns, a dictionary ===================
from collections import defaultdict
columns = defaultdict(list) #we want a list to append each value in each column to

with open(filename) as f:
    reader = csv.DictReader(f) #create a reader which represents rows in a dictionary form
    for row in reader: #this will read a row as {column1: value1, column2: value2,...}
        for (k,v) in row.items(): #go over each column name and value 
            columns[k].append(v) #append the value into the appropriate list based on column name k


# X is the input variable to the model regression
x1 = np.array(columns['Time'],dtype=np.float64)/60
del columns['Time']
x2 = np.array(columns.keys(),dtype=np.float64)
# ===============================================================================

# Main loop for phase 1: iterate over each column and perform regression
# Step 3. Parse a regression model based on some text string

# Create a parameters list to be used in estimation
params = Parameters()
params.add('k',   value= 0.05) # initial guess is 0.05 based on experience
params.add('A', value= 110) # initial guess is 110 based on previous fitting

# objective function to minimize, in this case this is minimizing residual
def phase1_residual(params, x, data):   
    return (phase1_model(params,x)- data)

def phase1_model(params, x):   
    k = params['k'].value
    A = params['A'].value
    
    model = A*(1-np.exp(-k*x)) # this is the model
    return model
        

result = dict()
# Main loop for phase 1 parameter estimation
i = 0 # counter for color purposes
P = dict()
pylab.close('all')


for (k,column_data) in columns.items():
    print "Now estimating the dataset %s" % k
    column_data = np.array(column_data ,dtype=np.float64)
    
    # Use nonlinear least squares regression to estimate the parameters
    result = minimize(phase1_residual, params, args=(x1, column_data))
    
    report_fit(params,show_correl=False)
    
    param_string = ""
    # store the estimated parameters in dictioary P for later use
    for (param_name,param_val) in params.items():
        param_string = param_string + "| %s: %0.4f| " % (param_name,param_val.value)
        if param_name in P.keys():
            P[param_name] = np.append(P[param_name],param_val.value)
        else:
            P[param_name]=np.array([param_val.value])

    
    final = column_data + result.residual
    # output the fit results in a figure
    try:       
        pylab.plot(x1, column_data, colors[i]+"o",label='Concentration (uM): '+k)
        x1_plot = np.linspace(0,max(x1),100)
        pylab.plot(x1_plot, phase1_model(params,x1_plot), colors[i],label='Fit:'+ param_string)
    except:
        pass
    i = i + 1

# label the graph
pylab.title('Phase 1 - Least Squares Regression Fit Results')
pylab.grid(True)
pylab.xlabel('Time (min)')
pylab.ylabel('Concentration (uM)')
pylab.legend()
pylab.show()

print ""
print """


    Phase 1 complete, the parameter estimates going into phase 2 are as follows...   
"""

for (k,v) in P.items():
    print k
    print v

print """


    Now starting phase 2 regression

    
"""
# Phase 2 - Assume the parameters themselves are also functions of 3 parameters

# objective function to minimize
def phase2_residual(params, x, P):
    return phase2_model(params,x)-P['k']


def phase2_model(params, x):
    IC50 = params['IC50'].value   
    v0 = params['v0'].value
    dv = params['dv'].value
    
    model = v0 - dv*(x/(x+IC50))   
    return model

p2_params = Parameters()
p2_params.add('IC50', value= 4)
p2_params.add('v0', value= 0.0075)
p2_params.add('dv', value= 0.0075)

phase2_results = minimize(phase2_residual, p2_params, args=(x2, P))
report_fit(p2_params,show_correl=False)

final = P['k'] + phase2_results.residual
p2_param_string = ""
P2 = dict()

# similar to phase 1, store results in dict and generate a string output
for (param_name,param_val) in p2_params.items():
    p2_param_string = p2_param_string + "|%s: %0.4f| " % (param_name,param_val.value)
    if param_name in P2.keys():
        P2[param_name] = np.append(P2[param_name],param_val.value)
    else:
        P2[param_name]=np.array([param_val.value])

try:
    pylab.figure()
    pylab.plot(x2, P['k'], "ro",label='Measured: '+k)
    x2_plot = np.linspace(0,20,50)
    pylab.plot(x2_plot, phase2_model(p2_params,x2_plot),"k", label='Fit: '+ p2_param_string)
except:
    pass

# formatin/labeling the output figure
pylab.grid(True)
pylab.xlabel('[dCTP] (uM)')
pylab.ylabel('Rate (min^-1)')
pylab.title('Phase 2 - Regression of hyperparameters based on phase 1 rate constants')
pylab.legend()
pylab.show()

