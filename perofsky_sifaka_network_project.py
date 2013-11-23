###Amanda Perofsky
###programming for biologists project
###Nov 2013
###note: this code is for Python 2.7
##note: some of percolation code was adapted from a teaching module by TJ Hladish for the Summer Institute in the Statistics and Modeling in Infectious Diseases

##A. determine the epidemic size distribution for the sifaka contact network

# Import functions
from networkx import *
from random import *
import pylab
import matplotlib
import matplotlib.pyplot as plt

# 1.) Define parameters
T = 0.5 # probability of transmission; this can be varied;would expect a high T value for an STD

# 2.) Construct network
G = networkx.read_edgelist('mating_edges.txt', nodetype = int) #read in edgelist txt file to create network object
print G.edges() #make sure all the edges were imported
net_size = float(G.order()) #return the number of nodes in the graph
results = [] #start off with empty results list

N, K = G.order(), G.size() #network size, number of edges
avg_deg = float(K)/N #average degree
print "Nodes: ", N
print "Edges: ", K
print "Average degree: ", avg_deg

###degree sequence plot
degree_sequence = sorted(networkx.degree(G).values(), reverse=True)
dmax= max(degree_sequence)

####other degree sequence code
deg_seq = G.degree().values()
print(deg_seq)

####plot G 
networkx.draw_spring(G)
plt.show()

# Average clustering coefficient
ccs = networkx.clustering(G)
avg_clust = sum(ccs.values()) / len(ccs)
print(avg_clust)

# 3.) Store node states
for i in range(500): #run the model 500 times 
    patient_zero = choice(G.nodes()) # Randomly choose one node to be patient zero
    infected = [patient_zero] #start off with patient zero as infected
    recovered = [] #start off with an empty list of recovered nodes
    
# 4.) Run the simulation
    while len(infected) > 0: #while the list of infected nodes is at least 1
        v = infected.pop(0) #v is the list of infected nodes (with the first item removed)
        for u in G.neighbors(v): #G.neighbors(v) returns a list of the nodes connected to the node v 
            if u not in infected and u not in recovered: #if node u is not infected or recovered (i.e., it is susceptible)
                if random() < T: #if a randomly generated number is less than the T value
                    infected.append(u) #infect u and append it to the infected nodes list
        recovered.append(v) #add v to the list of recovered nodes

    # 5.) Tally the total fraction that got infected
    epi_size = len(recovered)/net_size #epidemic size is length of recovered nodes list divided by network size
    results.append( epi_size ) #add epi_size outcome to the results list

plt.hist(results, bins = 10) #histogram of epidemic size
plt.xlabel("Total fraction of network infected")
plt.ylabel("Frequency")
plt.ylim([0,300])
plt.title("Distribution of epidemic sizes for 2008 mating season, T=0.5, 500 simulations")
plt.show() #waht does the distribution look like? poisson, scale-free,etc.?


####B. investigate the relationship between transmissability T and epidemic size S
###run percolation simulation many times for different values of T
###for each T value, average the outcomes (epidemic sizes)
###create a plot w/ the resulting mean epidemic sizes

net_size = float(G.order()) #return the number of nodes in the graph (network size)
results = [] #start off with an empty results list
mean_epi_sizes = [] ###start off with an empty epidemic size list

T_values =  [i/20. for i in range(21)] #generate a list of 20 T values from 0.0 to 1.0

####FIRST LOOP
for t in T_values:

    ####SECOND LOOP
    for i in range(100): #run the model 100 times 
        patient_zero = choice(G.nodes()) # Randomly choose one node to be patient zero
        infected = [patient_zero] #start off with patient zero as infected
        recovered = [] #start off with an empty list of recovered nodes
    
# 4.) Run the simulation
        while len(infected) > 0: #while the list of infected nodes is at least 1
            v = infected.pop(0) #v is the list of infected nodes (with the first item removed)
            for u in G.neighbors(v): #G.neighbors(v) returns a list of the nodes conntect to the node v 
            # We need to make sure Node u is still susceptible
                if u not in infected and u not in recovered: #if node u is not infected or recovered (i.e., it is susceptible)
                    if random() < t: #if a randomly generated number is less than the T value
                        infected.append(u) #infect u and append it to the infected nodes list
            recovered.append(v) #add v to the list of recovered nodes at the end of the simulation

    # 5.) Tally the total fraction that got infected
        epi_size = len(recovered)/net_size #epidemic size is length of recovered nodes list divided by network size (fraction of pop infected)
        results.append( epi_size ) #add epi_size outcome to the results list
        average_epi = float(sum(results))/len(results) #calculate the average epidemic size for the T value
        ####END SECOND LOOP

    ###FIRST LOOP
    mean_epi_sizes.append(average_epi) #create a list of average epidemic sizes from all the T values

print(mean_epi_sizes) #use this for plotting in R
print(T_values) ##use this for plotting in R 

from pylab import scatter,show
scatter(T_values, mean_epi_sizes) 
plt.xlabel("Transmissability")
plt.ylabel("Mean epidemic size")
plt.xlim([-0.1,1.1])
plt.title("Transmissability sensitivity analysis mating season, 100 simulations/T value")
show()

