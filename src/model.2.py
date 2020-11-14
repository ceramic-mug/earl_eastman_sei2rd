#########################################
# title: model.2                        #
# author: Joshua Eastman                #
# description: Dynamic SEIRD model      #
#     based on real-time geotemporal    #
#     params. Not your grandma's model. #
#########################################

## (1) Import Libraries

import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
import datetime

## (2) Import Flow Data
##     and compute populations

####################################################
# realFlowsMatrix.csv                              
#
#   from        to          n
#   a           a           x
#   a           b           x
#   a           c           x
#   ...         ...         ...
#
# number of rows (no header) = (number of nodes)^2
#   "from" and "to" ordering must be identical
####################################################
flow_data = pd.read_csv('../in/realFlowsMatrix.csv')

# Extract a list of node names based on the dimensions of flow_data
node_names = flow_data['to'].values[0:int(math.sqrt(len(flow_data['to'].values)))]

# Create a matrix using flow_data n values
# Note: F required for filling in from-to order
flow_matrix = flow_data['n'].values
flow_matrix = np.reshape(flow_matrix, (len(node_names),len(node_names)), order='F')

# Extract the population of each node
pop = [sum(flow_matrix[i,:]) for i in range(len(flow_matrix[:,1]))]

# total population
pop_total = sum(pop)

# fractional population
pop_frac = pop/pop_total

# largest and smallest nodes (indexes)
index_largest_node = np.argmax(pop)
index_smallest_node = np.argmin(pop)

## (3) Set and Solve Biological Parameters

R0 = 2.5 # number of secondary infections per
         # primary infection in a totally
         # naive population

mean_infectious_period = 5 # days
gamma = 1 / mean_infectious_period

mean_latent_period = 3 # days
sigma = 1 / mean_latent_period

proportion_symptomatic = 0.85
omega = proportion_symptomatic

proportion_severe = 0.05
delta = proportion_severe

natural_lifespan = 70*365 # days
mu = 1/natural_lifespan # 1/day

beta = R0 * (gamma + delta + mu)

## (4) Construct Data Structure

# time array
dt = 1 # day
time = np.arange(start=0, stop=365, step = dt)

# model categories and subpopulations
subpopulations = ['home','commuter']
categories = ['S','E','Is','Ia','R','D']

# data structure to hold model output
# dimensions: #nodes x timesteps x 2 x 6
model_out = np.empty(shape=(len(node_names),len(time),len(subpopulations),len(categories)))

## (5) Population Partitioning and Beta Construction

beta_arr = np.empty(shape=(len(pop),2))
for i in range(len(pop)):
    home = flow_matrix[i,i]
    comm = pop[i]-home
    model_out[i,0,0,0] = home # diagonal
    model_out[i,0,1,0] = comm # row - diagonal

    # calculate betas for home and commuter classes
    beta_arr[i,0] = beta / (comm + sum(flow_matrix[:,i]))

    comm_divisor = comm
    for j in range(len(pop)):
        if flow_matrix[i,j] != 0:
            comm_divisor += flow_matrix[j,j] * (flow_matrix[i,j] / comm)

    beta_arr[i,1] = beta /  comm_divisor

# set all other categories equal to 0 at start of run

for i in range(len(categories)-1):
    model_out[:,:,:,i+1] = 0

## (6) Functions for constructing transition matrices

# t is previous timestep, n is node
def I_mixing(t, n):
    I_h = 0
    I_c = 0
    for m in range(len(pop)):
        if m != n:
            I_h += flow_matrix[m,n]/model_out[m,0,1,0] * sum([model_out[m,t,1,2],model_out[m,t,1,3]])
            I_c += flow_matrix[n,m]/model_out[n,0,1,0] * sum([model_out[m,t,0,2],model_out[m,t,0,3]])
        else:
            I_h += sum([model_out[m,t,0,2],model_out[m,t,0,3]])
            I_c += sum([model_out[m,t,0,2],model_out[m,t,0,3]])
    return(I_h, I_c)

def transition(t,n):
    I_set = I_mixing(t-1,n)
    
    h_S = np.array([-(beta_arr[n,0] * I_set[0]), mu, mu, mu, mu, 0])
    h_E = np.array([(beta_arr[n,0] * I_set[0]), - (sigma + mu), 0, 0, 0, 0])
    h_Is = np.array([0, omega * sigma, -(mu + (delta * gamma) + gamma), 0, 0, 0])
    h_Ia = np.array([0, (1-omega) * sigma, 0, -(mu + gamma), 0, 0])
    h_R = np.array([0, 0, gamma, gamma, -mu, 0])
    h_D = np.array([0, 0, delta * gamma, 0, 0, 0])

    h_matrix = np.matrix([h_S, h_E, h_Is, h_Ia, h_R, h_D])

    c_S = np.array([-(beta_arr[n,1] * I_set[1]), mu, mu, mu, mu, 0])
    c_E = np.array([(beta_arr[n,1] * I_set[1]), - (sigma + mu), 0, 0, 0, 0])
    c_Is = np.array([0, omega * sigma, -(mu + (delta * gamma) + gamma), 0, 0, 0])
    c_Ia = np.array([0, (1-omega) * sigma, 0, -(mu + gamma), 0, 0])
    c_R = np.array([0, 0, gamma, gamma, -mu, 0])
    c_D = np.array([0, 0, delta * gamma, 0, 0, 0])

    c_matrix = np.matrix([c_S, c_E, c_Is, c_Ia, c_R, c_D])

    return (h_matrix, c_matrix)


## (7) Seed the Model with infected 

model_out[index_largest_node,0,0,2] = 1
model_out[index_largest_node,0,0,0] -= 1

## (8) MAIN LOOP

for timestep in time:
    if timestep == 0:
        continue

    for n in range(len(pop)):

        t_matrices = transition(timestep, n)

        # home
        model_out[n, timestep, 0, :] =  model_out[n, timestep-1, 0, :] + dt * t_matrices[0].dot(model_out[n, timestep-1, 0, :])
        # comm
        model_out[n, timestep, 1, :] =  model_out[n, timestep-1, 1, :] + dt * t_matrices[1].dot(model_out[n, timestep-1, 1, :])

## (9) Save to CSV

node_out = [None]*len(pop)*len(time)
time_out = [None]*len(pop)*len(time)
S_out = [None]*len(pop)*len(time)
E_out = [None]*len(pop)*len(time)
Is_out = [None]*len(pop)*len(time)
Ia_out = [None]*len(pop)*len(time)
R_out = [None]*len(pop)*len(time)
D_out = [None]*len(pop)*len(time)

for n in range(len(pop)):
    for timestep in time:
        node_out[n * len(time) + timestep] = node_names[n]
        time_out[n * len(time) + timestep] = timestep
        S_out[n * len(time) + timestep] = sum(model_out[n,timestep,:,0])
        E_out[n * len(time) + timestep] = sum(model_out[n,timestep,:,1])
        Is_out[n * len(time) + timestep] = sum(model_out[n,timestep,:,2])
        Ia_out[n * len(time) + timestep] = sum(model_out[n,timestep,:,3])
        R_out[n * len(time) + timestep] = sum(model_out[n,timestep,:,4])
        D_out[n * len(time) + timestep] = sum(model_out[n,timestep,:,5])

data_dict = {'node':node_out,
             't':time_out,
             'S':S_out,
             'E':E_out,
             'Is':Is_out,
             'Ia':Ia_out,
             'R':R_out,
             'D':D_out}

output_datafile = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
data_out = pd.DataFrame(data = data_dict)
data_out.to_csv('../out/'+output_datafile+'.csv',index=False)