import networkx as nx
import numpy as np
import itertools
from scipy.spatial import distance
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import csv
import pdb
import globals as const
from funcs import *
from math import isclose

# Constants for simulation
dt = const.dt
#var_dt = True

# dimensions of the cell 
l_apical = const.l_apical 

# Set the arcs list
inner_arc = const.inner_arc
outer_arc = const.outer_arc

# mechanical parameters
l0_apical = l_apical
mu_apical = const.mu_apical         
myo_beta = const.myo_beta 
eta = const.eta 
l_mvmt = const.l_mvmt

# initialize the tissue
G, centers, num_api_nodes, circum_sorted, belt, triangles = tissue_2d()
pit_centers = const.pit_centers 

# Starting from t=0
t = 0
num_inter = 0 
blacklist = [] 
contract = [True for counter in range(0,num_inter)]

# t=initial nx Graph in pickled form for plotting later
print(t) 
file_name = 't' + str(int(t)) 
nx.write_gpickle(G,file_name + '.pickle')
np.save(file_name,circum_sorted) 

while t <= const.t_final:
    
    # increment t by dt
    # initialize force_dict back to zeros
    t = round(t+dt,1)
    print(dt, t) 
    pos = nx.get_node_attributes(G,'pos')
    force_dict = {new_list: np.zeros(2,dtype=float) for new_list in G.nodes()} 

    # Update myosin on a fictitious pit (no resemblance to SG geometry)
    if t == const.t_pit: 
        for node in pit_centers: 
            if node == 0:
                myo = 1.39*const.pit_strength
            for neighbor in G.neighbors(node): 
                G[node][neighbor]['myosin'] = const.pit_strength 
        print("Pit is established")

#    if t > const.t_intercalate:
#        if contract[0] == True:
#            G[301][302]['myosin'] = const.belt_strength*(t-const.t_intercalate) 
    
    # update myosin on inner arc 
#    if t == const.t_1:
#        for i in range(0,len(inner_arc)):
#            G[inner_arc[i-1]][inner_arc[i]]['myosin'] = const.belt_strength     
#        print("Inner arc established")

#    # update myosin on outer arc 
#    if t == const.t_2:
#        for i in range(0,len(outer_arc)):
#            G[outer_arc[i-1]][outer_arc[i]]['myosin'] = const.belt_strength     
#        print("Outer arc established")

#    # update myosin on belt
#    if t == const.t_belt:
#        for i in range(0,len(belt)):
#            G[belt[i-1]][belt[i]]['myosin'] = const.belt_strength     
#        print("Belt established") 

    for node in G.nodes(): 
        # update force on each node  
        force = [0.0,0.0]
    
        # Elastic forces due to the cytoskeleton 
        for neighbor in G.neighbors(node):
            a = pos[node]
            b = pos[neighbor]
            
            dist = distance.euclidean(a,b)
            direction = unit_vector(a,b)
            
            if (node in inner_arc) and (neighbor in inner_arc):
                magnitude = elastic_force(dist, G[node][neighbor]['l_rest'], const.mu_apical) 
            else:
                magnitude = elastic_force(dist, G[node][neighbor]['l_rest'], mu_apical) 
            force = np.sum([force,magnitude*np.array(direction)],axis=0)
            
            # Force due to myosin
            magnitude = myo_beta*G[node][neighbor]['myosin']
            force = np.sum([force, magnitude*np.array(direction)],axis=0)

        force_dict[node] = np.add(force_dict[node], force) 
   
    # update location of node 
    pos = nx.get_node_attributes(G,'pos')
    
    for node in force_dict:
        G.node[node]['pos'] = d_pos(pos[node],force_dict[node],dt)

# Save nx Graph in pickled form for plotting later
    
    if t % 1 == 0: 
        file_name = 't' + str(round(t)) 
        nx.write_gpickle(G,file_name + '.pickle')
        np.save(file_name,circum_sorted)
