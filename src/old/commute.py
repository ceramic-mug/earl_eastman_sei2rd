# to generate numbers of people who travel from each NYC
# borough to every other NYC borough and people who travel
# out of NYC

# import packages
import pandas as pd
import numpy as np

# read in our data
c2015 = pd.read_excel('/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/dat/commute_2015.xlsx')
pertinent = pd.read_csv('/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/dat/commute_counties.csv')

cols = ['from','to','n']

# dictionary to hold our commuter computed values
commute = pd.DataFrame(columns=cols)

for index, row in pertinent.iterrows():
    in_n = 0
    out_n = 0
    from_dat = c2015.query('from_state==\"'+row['state']+'\" & from_county==\"'+row['county']+'\" ')
    for i, r in pertinent.iterrows():
        in_n += from_dat[from_dat.to_county==r['county']].n.tolist()[0]
        a = pd.DataFrame(data=[[row['county'],r['county'],from_dat[from_dat.to_county==r['county']].n.tolist()[0]]],columns=cols)
        commute = commute.append(a)
    out_n = sum(from_dat['n']) - in_n
    commute = commute.append(pd.DataFrame(data=[[row['county'],'out',out_n]],columns=cols))

print(commute.head())

commute.to_csv('/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/dat/pertinent_commute.csv',index=False)
