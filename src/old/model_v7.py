import numpy as np
import pandas as pd
import re
import random
import datetime
import os

##### Useful regular expressions

color_getter = re.compile('#.*')

###### Reference lists and dictionaries ######

categories = [
'm',
'c',
'i',
'h'
]

boroughs = [
'Bronx',
'Manhattan',
'Brooklyn',
'Staten Island',
'Queens'
]

county_to_borough = {
    'Bronx County': 'Bronx',
    'New York County': 'Manhattan',
    'Kings County': 'Brooklyn',
    'Richmond County': 'Staten Island',
    'Queens County': 'Queens',
    'out':'out'
}

###### functions ######

# convert feet to kilometers
def ft_to_m(_d_):
    return 0.3048 * _d_

###### import data ######

# import commuter flows data
commuting_flows_by_borough = pd.read_csv('/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/dat/pertinent_commute.csv')
# change names to common names
commuting_flows_by_borough['from'] = [county_to_borough[i] for i in commuting_flows_by_borough['from']]
commuting_flows_by_borough['to'] = [county_to_borough[i] for i in commuting_flows_by_borough['to']]

# import nta population data
nta_populations = pd.read_csv('/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/dat/New_York_City_Population_By_Neighborhood_Tabulation_Areas.csv')
nta_populations.head()
nta_populations = nta_populations[nta_populations['Year']==2010].sort_values(by=['NTA Code'])

# import nta distances to nearest subway stop
nta_distances_to_subway = pd.read_csv('/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/dat/nta_subway_stop_distances_to_closest.csv')
# sort by NTA index
nta_distances_to_subway = nta_distances_to_subway.sort_values(by=['nta_code'])

# import subways per Borough
borough_subways = pd.read_csv('/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/dat/borough_subways.csv')
subways = np.unique(borough_subways['color'].values)
# get a list of non-zero population NTAs for iteration
ntas = [i for i in nta_distances_to_subway['nta_code'].values.tolist() if nta_populations[nta_populations['NTA Code']==i]['Population'].values[0] != 0]

##### END IMPORT DATA #####

##### Compute Paramaters ######

#### PARAM 1: NTA sub_boxes populations

# sum population data by borough
borough_populations = {}

# use NTA populations to create borough populations for proportional commuter flow
for borough in nta_populations.Borough.unique():
    borough_populations[borough] = sum([i for i in nta_populations[nta_populations['Borough']==borough]['Population']])

# of those who commute, what proportion commutes where?
commuting_flows_by_borough['proportion'] = [row['n']/sum([i for i in commuting_flows_by_borough[commuting_flows_by_borough['from']==row['from']]['n']]) for index, row in commuting_flows_by_borough.iterrows()]

# from research, what percentage of the population commutes/works?
prop_commute = 0.648 # 1NYC and vital statistics of NY state, 2007

# then the number who stays at home/doesn't work is
prop_home = 1 - prop_commute

# create a dictionary to hold the sums of proportion of people who commute into other boroughs per borough
inter_borough_commuter_proportions = {}
for borough in boroughs:
    inter_borough_commuter_proportions[borough] = sum([i for i in commuting_flows_by_borough.loc[((commuting_flows_by_borough['from']==borough) & (commuting_flows_by_borough['to']!='out') & (commuting_flows_by_borough['to']!=borough))]['proportion']])

# add a column to the subway distances dataset that is the distance in kilometers
nta_distances_to_subway['distance(m)'] = [ft_to_m(i) for i in nta_distances_to_subway['distance(ft)']]

# function for subway ridership per population (_dist_ in meters, _borough_ string)
def subway_prop_population(_dist_,_borough_):

    a =  inter_borough_commuter_proportions[_borough_] * prop_commute * (-0.00006 * _dist_ + 1) # up from 0.53 at d=0 from Gutierrez et al, 2011

    # can't have a negative proportion
    if a <= 0:
        a = 0

    return a

# function for car ridership per population (_dist_ in meters)
def car_prop_population(_dist_,_borough_):
    return inter_borough_commuter_proportions[_borough_] * prop_commute - subway_prop_population(_dist_,_borough_)

# function for those who work within the borough
def in_commute_prop_population(_dist_,_borough_):
    return prop_commute * commuting_flows_by_borough.loc[((commuting_flows_by_borough['from']==_borough_) & (commuting_flows_by_borough['to']==_borough_))]['proportion'].values[0]

# Each NTA will have 4 sub-boxes:
#   - inter-borough commuters by subway
#   - inter-borough commuters by car
#   - home borough commuters
#   - stay-at-home or commute out of city
# The populations for each of these boxes will be computed for each NTA by the following function:

def partition(_dist_,_borough_,_pop_):
    m_ = int(subway_prop_population(_dist_,_borough_)*_pop_)
    c_ = int(car_prop_population(_dist_,_borough_)*_pop_)
    i_ = int(in_commute_prop_population(_dist_,_borough_)*_pop_)
    h_ = int(_pop_ - (m_+c_+i_))
    return [m_,c_,i_,h_]

#### PARAM 2: Computing effective population coefficients for subway lines, NTAs, and Boroughs

# find the proportion of each metro, car box goes to the borough in question (first nesting) from each other borough (second nesting)
borough_propto = {}
for _borough_ in boroughs:
    borough_propto[_borough_] = {}
    _other_boroughs_ = [i for i in boroughs if i!=_borough_]
    for b in _other_boroughs_:
        borough_propto[_borough_][b] = commuting_flows_by_borough.loc[((commuting_flows_by_borough['from']==b)&(commuting_flows_by_borough['to']==_borough_))]['proportion'].values[0]/inter_borough_commuter_proportions[b]


# compute the effective borough I/N at each timestep
def effective_borough(_borough_,_nta_dict_,_model_,_timestep_):
    # which NTAs send people into this borough
    _outer_ntas_ = [i for i in nta_populations.loc[nta_populations['Borough']!=_borough_]['NTA Code'].values if i in ntas]
    _n_ = 0
    _i_ = 0

    # add from people commuting in from the outside
    for _nta_ in _outer_ntas_:
        if sum(sum(_model_[_timestep_,_nta_dict_[_nta_],:2,:-1])) != np.nan and sum(sum(_model_[_timestep_,_nta_dict_[_nta_],:2,:-1])) != 0:
            _n_ += borough_propto[_borough_][nta_populations[nta_populations['NTA Code']==_nta_]['Borough'].values[0]]*sum(sum(_model_[_timestep_,_nta_dict_[_nta_],:2,:-1]))
            _i_ += borough_propto[_borough_][nta_populations[nta_populations['NTA Code']==_nta_]['Borough'].values[0]]*sum(sum(_model_[_timestep_,_nta_dict_[_nta_],:2,2:4]))

    _inner_ntas_ = [i for i in nta_populations.loc[nta_populations['Borough']==_borough_]['NTA Code'].values if i in ntas]

    # add from people commuting internally
    for _nta in _inner_ntas_:
        if sum(_model_[_timestep_,_nta_dict_[_nta],2,:-1]) != np.nan and sum(_model_[_timestep_,_nta_dict_[_nta],2,:-1]) != 0:
            _n_ += sum(_model_[_timestep_,_nta_dict_[_nta],2,:-1])
            _i_ += sum(_model_[_timestep_,_nta_dict_[_nta],2,2:4])
    return (_i_,_n_)

# compute the effective NTA I/N at each timestep
def effective_nta(_nta_,_nta_dict_,_model_,_timestep_):
    _i_ = sum(sum(_model_[_timestep_,_nta_dict_[_nta_],:,2:4]))
    _n_ = sum(sum(_model_[_timestep_,_nta_dict_[_nta_],:,:-1]))
    return (_i_,_n_)

# compute the number of distinct subway lines (colors) an NTA subway rider population rides
def nta_colors(_nta_code_):
    _colors_ = [i for i in nta_distances_to_subway[nta_distances_to_subway['nta_code']==_nta_code_][['color1','color2','color3','color4']].values[0].tolist() if not pd.isnull(i)]
    _unique_ = [i for i in np.unique(np.array(_colors_))]
    _n_ = len(_unique_)
    return (_n_, _unique_)

# compute the effective train line I/N at each timestep
def effective_train(_train_color_,_nta_dict_,_model_,_timestep_):
    _n_ = 0
    _i_ = 0

    for _nta_ in ntas:
        a = nta_colors(_nta_)

        # if the nta's subway population rides this subways
        if _train_color_ in a[1]:
            _n_ += sum(_model_[_timestep_,_nta_dict_[_nta_],0,:-1])/a[0]
            _i_ += sum(_model_[_timestep_,_nta_dict_[_nta_],0,2:4])/a[0]

        # otherwise, use the commuter influxes
        else:
            for _borough_ in borough_subways[borough_subways['color']==_train_color_]['borough'].values:
                if _nta_ not in nta_populations[nta_populations['Borough']==_borough_]['NTA Code'].values:
                    _n_ += borough_propto[_borough_][nta_populations[nta_populations['NTA Code']==_nta_]['Borough'].values[0]]*sum(_model_[_timestep_,_nta_dict_[_nta_],0,:-1])/len(borough_subways[borough_subways['borough']==_borough_].values)
                    _i_ += borough_propto[_borough_][nta_populations[nta_populations['NTA Code']==_nta_]['Borough'].values[0]]*sum(_model_[_timestep_,_nta_dict_[_nta_],0,2:4])/len(borough_subways[borough_subways['borough']==_borough_].values)

    return (_i_,_n_)

######################## END PARAMS

##### MATRICES #######

mean_latent_period = 3
proportion_symptomatic = 0.86834
mean_infectious_period = 5
proportion_severe_cases = 0.05

params = {'mean latent period': mean_latent_period,
          'proportion symptomatic': proportion_symptomatic,
          'mean infectious period': mean_infectious_period,
          'proportion severe cases': proportion_severe_cases
          }

# effective populations a commuter from a given NTA is exposed to
def commuter_eff(_nta_,_effective_dict_):

    term = 0
    home_borough = nta_populations[nta_populations['NTA Code']==_nta_]['Borough'].values[0]
    # ticker
    i = 0
    # total encountered infected population
    _i_ = 0
    # total population
    _n_ = 0

    for _borough_ in borough_propto[home_borough].keys():
        # ticker ticks up one
        i += 1
        # add to the total infected population the computed effective infected population of the borough
        _i_ += borough_propto[home_borough][_borough_] * _effective_dict_['borough'][_borough_][0]
        # add to the total population the effective population of the borough
        _n_ += borough_propto[home_borough][_borough_] * _effective_dict_['borough'][_borough_][1]

    return (_i_,_n_)

def metro_eff(_nta_,_effective_dict_):

    home_borough = nta_populations[nta_populations['NTA Code']==_nta_]['Borough'].values[0]
    # initialize effective infected and total populations encountered
    _i_ = 0
    _n_ = 0

    for _borough_ in borough_propto[home_borough].keys():

        # averages for the number of lines in each borough and the proportion of the commuting population going to that borough

        _i_ += borough_propto[home_borough][_borough_] * np.average(np.array([_effective_dict_['metro'][i][0] for i in borough_subways[borough_subways['borough']==_borough_]['color'].values]))
        _n_ += borough_propto[home_borough][_borough_] * np.average(np.array([_effective_dict_['metro'][i][1] for i in borough_subways[borough_subways['borough']==_borough_]['color'].values]))

    return (_i_,_n_)

def home_eff(_nta_,_effective_dict_):
    return (_effective_dict_['nta'][_nta_][0],_effective_dict_['nta'][_nta_][1])

def inborough_eff(_nta_,_effective_dict_):
    return (_effective_dict_['borough'][nta_populations[nta_populations['NTA Code']==_nta_]['Borough'].values[0]][0],_effective_dict_['borough'][nta_populations[nta_populations['NTA Code']==_nta_]['Borough'].values[0]][1])

# standard beta and i and n
def home_term(beta_nta_,_nta_,_effective_dict_):
    return beta_nta_ * home_eff(_nta_,_effective_dict_)[0]/home_eff(_nta_,_effective_dict_)[1]

# weighted proportional sum of the beta values for home and for borough multiplied by the effective i over n total that we're exposed to
def inborough_commuter_term(beta_nta_,beta_borough_,_nta_,_effective_dict_):
    _i_ = inborough_eff(_nta_,_effective_dict_)[0] + home_eff(_nta_,_effective_dict_)[0]
    _n_ = inborough_eff(_nta_,_effective_dict_)[1] + home_eff(_nta_,_effective_dict_)[1]
    prop_borough = inborough_eff(_nta_,_effective_dict_)[1] / _n_
    prop_home = home_eff(_nta_,_effective_dict_)[1] / _n_
    weighted_beta = sum([prop_home*beta_nta_,prop_borough*beta_borough_])
    return weighted_beta * _i_ / _n_

# weighted proportional sum of the beta values multiplied by the populations exposed to
def car_commuter_term(beta_nta_,beta_borough_,_nta_,_effective_dict_):
    _i_ = commuter_eff(_nta_,_effective_dict_)[0] + home_eff(_nta_,_effective_dict_)[0]
    _n_ = commuter_eff(_nta_,_effective_dict_)[1] + home_eff(_nta_,_effective_dict_)[1]
    prop_borough = commuter_eff(_nta_,_effective_dict_)[1] / _n_
    prop_home = home_eff(_nta_,_effective_dict_)[1] / _n_
    weighted_beta = sum([prop_home*beta_nta_,prop_borough*beta_borough_])
    return weighted_beta * _i_ / _n_

# weighted proportional sum of the beta values multiplied by the populations exposed to
def metro_commuter_term(beta_nta_,beta_metro_,beta_borough_,_nta_,_effective_dict_):
    _i_ = metro_eff(_nta_,_effective_dict_)[0] + commuter_eff(_nta_,_effective_dict_)[0] + home_eff(_nta_,_effective_dict_)[0]
    _n_ = metro_eff(_nta_,_effective_dict_)[1] + commuter_eff(_nta_,_effective_dict_)[1] + home_eff(_nta_,_effective_dict_)[1]
    prop_metro = metro_eff(_nta_,_effective_dict_)[1] /_n_
    prop_borough = commuter_eff(_nta_,_effective_dict_)[1] / _n_
    prop_home = home_eff(_nta_,_effective_dict_)[1] / _n_
    weighted_beta = sum([prop_home*beta_nta_,prop_borough*beta_borough_,prop_metro*beta_metro_])
    return weighted_beta * _i_ / _n_

# defined according to array (S,E,Is,Ia,R,D)
def matrices(beta_metro_,beta_nta_,beta_borough_,_effective_dict_,_nta_):

    home_matrix = np.array([[-(home_term(beta_nta_,_nta_,_effective_dict_)),0,0,0,0,0],
                            [home_term(beta_nta_,_nta_,_effective_dict_),-mean_latent_period**(-1),0,0,0,0],
                            [0,proportion_symptomatic*mean_latent_period**(-1),-mean_infectious_period**(-1)*(proportion_severe_cases+1),0,0,0],
                            [0,(1-proportion_symptomatic)*mean_latent_period**(-1),0,-mean_infectious_period**(-1),0,0],
                            [0,0,mean_infectious_period**(-1),mean_infectious_period**(-1),0,0],
                            [0,0,mean_infectious_period**(-1)*proportion_severe_cases,0,0,0]])

    subway_matrix = np.array([[-(metro_commuter_term(beta_nta_,beta_metro_,beta_borough_,_nta_,_effective_dict_)),0,0,0,0,0],
                            [metro_commuter_term(beta_nta_,beta_metro_,beta_borough_,_nta_,_effective_dict_),-mean_latent_period**(-1),0,0,0,0],
                            [0,proportion_symptomatic*mean_latent_period**(-1),-mean_infectious_period**(-1)*(proportion_severe_cases+1),0,0,0],
                            [0,(1-proportion_symptomatic)*mean_latent_period**(-1),0,-mean_infectious_period**(-1),0,0],
                            [0,0,mean_infectious_period**(-1),mean_infectious_period**(-1),0,0],
                            [0,0,mean_infectious_period**(-1)*proportion_severe_cases,0,0,0]])

    car_matrix = np.array([[-(car_commuter_term(beta_nta_,beta_borough_,_nta_,_effective_dict_)),0,0,0,0,0],
                            [car_commuter_term(beta_nta_,beta_borough_,_nta_,_effective_dict_),-mean_latent_period**(-1),0,0,0,0],
                            [0,proportion_symptomatic*mean_latent_period**(-1),-mean_infectious_period**(-1)*(proportion_severe_cases+1),0,0,0],
                            [0,(1-proportion_symptomatic)*mean_latent_period**(-1),0,-mean_infectious_period**(-1),0,0],
                            [0,0,mean_infectious_period**(-1),mean_infectious_period**(-1),0,0],
                            [0,0,mean_infectious_period**(-1)*proportion_severe_cases,0,0,0]])

    inborough_matrix = np.array([[-(inborough_commuter_term(beta_nta_,beta_borough_,_nta_,_effective_dict_)),0,0,0,0,0],
                            [inborough_commuter_term(beta_nta_,beta_borough_,_nta_,_effective_dict_),-mean_latent_period**(-1),0,0,0,0],
                            [0,proportion_symptomatic*mean_latent_period**(-1),-mean_infectious_period**(-1)*(proportion_severe_cases+1),0,0,0],
                            [0,(1-proportion_symptomatic)*mean_latent_period**(-1),0,-mean_infectious_period**(-1),0,0],
                            [0,0,mean_infectious_period**(-1),mean_infectious_period**(-1),0,0],
                            [0,0,mean_infectious_period**(-1)*proportion_severe_cases,0,0,0]])

    return {'m':subway_matrix,'c':car_matrix,'i':inborough_matrix,'h':home_matrix}

# compute the effectvie I/N coefficients for each thing in each category
def compute_effectives(_nta_dict_,_model_,_timestep_):
    effective_dict = {'borough':{},'nta':{},'metro':{}}

    for _borough_ in boroughs:
        effective_dict['borough'][_borough_] = effective_borough(_borough_,_nta_dict_,_model_,_timestep_)
    for _nta_ in ntas:
        effective_dict['nta'][_nta_] = effective_nta(_nta_,_nta_dict_,_model_,_timestep_)
    for _line_ in subways:
        effective_dict['metro'][_line_] = effective_train(_line_,_nta_dict_,_model_,_timestep_)

    return effective_dict

def main(time,initial_beta_metro,initial_beta_nta,initial_beta_borough,lockdown_factor):

    # initialize a 4-D array
    # dim
    # 1: time
    # 2: NTA
    # 3: category
    # 4: SIR box
    model = np.zeros([time,len(ntas),4,6])

    # map the NTA names to numbers in array
    nta_dict = {}

    # for all the NTAs with nonzero populations
    for i in range(len(ntas)):

        # the name is assigned to the number
        nta_dict[ntas[i]] = i

        # partition the NTA into its metro, car, in-borough, and home boxes
        vals = partition(nta_distances_to_subway[nta_distances_to_subway['nta_code']==ntas[i]]['distance(m)'].values[0],nta_populations[nta_populations['NTA Code']==ntas[i]]['Borough'].values[0],nta_populations[nta_populations['NTA Code']==ntas[i]]['Population'].values[0])

        # assign those initial populations to the susceptible box
        for j in range(len(vals)):
            model[0,i,j,0] = vals[j]

    nonzero_metro_ntas = [i for i in ntas if model[0,nta_dict[i],0,0] != 0]
    # seed the model with 4 random infected persons
    par = open("/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/out/notes-"+datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")+'.txt',"w+")
    par.write('Model paramaters:')
    for i in range(4):
        starter_nta = random.randint(0,len(nonzero_metro_ntas))
        model[0,nta_dict[nonzero_metro_ntas[starter_nta]],0,0] -= 1
        model[0,nta_dict[nonzero_metro_ntas[starter_nta]],0,3] += 1
        par.write('\nStarter NTA '+str(i)+' is '+nonzero_metro_ntas[starter_nta])

    # set the initial beta_nta value
    beta_metro = initial_beta_metro
    beta_nta = initial_beta_nta
    beta_borough = initial_beta_borough
    # write a file containing the paramaters
    params['initial beta metro'] = beta_metro
    params['initial beta nta'] = beta_nta
    params['initial beta borough'] = beta_borough
    for param in params.keys():
        par.write('\n'+param+': '+str(params[param]))
    par.write('\nLockdown factor: '+str(lockdown_factor))
    par.close()

    output_datafile = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    output_path = "/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/out/" + output_datafile +'.csv'
    # data output file
    f = open(output_path,"w+")
    f.write(','.join(['step','nta','cat','s','e','is','ia','r','d']))

    # write initial timestep with susceptible populations
    for nta in ntas:
        for i in range(len(categories)):
            f.write('\n'+','.join([str(0),nta,categories[i],','.join([str(i) for i in model[0,nta_dict[nta],i,:]])]))

    # timestep through the algorithm
    for step in range(time-1):

####### GOVERNMENT INTERVENTION (BASED ON NYC REALITY) ##########

#### two week lag

        # 20 percent metro reduction on day 7 (stay home)
        if step==20:
            for nta in ntas:
                not_riding_anymore = 0.2 * model[step,nta_dict[nta],0,:-1]
                # people leave the metro boxes
                model[step,nta_dict[nta],0,:-1] = model[step,nta_dict[nta],0,:-1] - not_riding_anymore
                # and add to the car boxes
                model[step,nta_dict[nta],1,:-1] = model[step,nta_dict[nta],1,:-1] + not_riding_anymore

        # 75% commuter population reduction on day 20 (stay at home order) --> stay home
        if step==33:
            for nta in ntas:
                for cat in range(len(categories)-1):
                    staying_home = 0.75 * model[step,nta_dict[nta],cat,:-1]

                    # leave their commuter box
                    model[step,nta_dict[nta],cat,:-1] = model[step,nta_dict[nta],cat,:-1] - staying_home
                    # add to the home box
                    model[step,nta_dict[nta],3,:-1] = model[step,nta_dict[nta],3,:-1] + staying_home

            # All the betas are reduced by a factor of five due to social distancing
            beta_nta = lockdown_factor*beta_nta
            beta_borough = lockdown_factor*beta_borough
            beta_metro = lockdown_factor*beta_metro

        # after day 7, symptomatic people stay home
        if step >= 20:
            for nta in ntas:
                for cat in range(len(categories)-1):
                    staying_home = model[step,nta_dict[nta],cat,2]
                    model[step,nta_dict[nta],cat,2] -= staying_home
                    model[step,nta_dict[nta],3,2] += staying_home

###################################################################

        # compute the effective borough, train line, and NTC coefficients
        #### CURRENTLY ALLOWS FOR METRO BETA INPUT #######
        effectives = compute_effectives(nta_dict,model,step)
        for nta in ntas:
            # create their matrix
            m = matrices(beta_metro,beta_nta,beta_borough,effectives,nta)
            for i in range(len(categories)):
                print(','.join([str(step),nta,categories[i],','.join([str(int(round(i))) for i in model[step,nta_dict[nta],i,:]])]))

                # compute
                model[step+1,nta_dict[nta],i,:] = [j for j in model[step,nta_dict[nta],i,:] + m[categories[i]].dot(model[step,nta_dict[nta],i,:])]

                f.write('\n'+','.join([str(step+1),nta,categories[i],','.join([str(i) for i in model[step+1,nta_dict[nta],i,:]])]))

    f.close()
    return output_datafile
