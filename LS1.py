import os
from sys import argv
import math
import time
import random


#######################################
#the class object for Simulated Annealing
#######################################
class SimAnneal(object):
    def __init__(self, coords, T=-1, alpha=-1, stopping_T=-1, stopping_iter=-1,restarts=-1,format_type=-1,random_seed=-1,time_cutoff=-1):
        self.coords = coords
        self.N = len(coords)
        self.format_type='EUC_2D' if format_type == -1 else format_type

        if random_seed == -1:
            random.seed()
        else:
            random.seed(random_seed)

        #self.T = math.sqrt(self.N) if T == -1 else T
        #self.stopping_temperature = 0.00000001 if stopping_T == -1 else stopping_T

        self.begin_T = 1e+10 if T == -1 else T
        self.T = 1e+10 if T == -1 else T
        self.stopping_temperature = 0.0001 if stopping_T == -1 else stopping_T

        #self.alpha = 0.995 if alpha == -1 else alpha
        self.alpha = 0.999 if alpha == -1 else alpha

        self.stopping_iter = 10000000 if stopping_iter == -1 else stopping_iter
        #self.stopping_iter = 100000 if stopping_iter == -1 else stopping_iter
        self.iteration = 1

        self.restarts= 10000 if restarts == -1 else restarts      # random restarts
        #self.restarts= 1 if restarts == -1 else restarts      # random restarts

        self.dist_matrix = self.to_dist_matrix(coords)
        self.nodes = [i for i in range(self.N)]

        self.cur_solution = self.initial_solution()        # Initialize greedy solution 
        #self.cur_solution = list(self.nodes)

        self.best_solution = list(self.cur_solution)

        self.cur_fitness = self.fitness(self.cur_solution)
        self.initial_fitness = self.cur_fitness
        self.best_fitness = self.cur_fitness

        self.fitness_list = [self.cur_fitness]

        self.best_iteration_solution=[]
        self.best_iteration_fitness=[]

        self.start_time=0 #initialize
        self.time_cutoff=20 if time_cutoff == -1 else time_cutoff #initialize
        self.time_fitness_list=[] #initialize

    def initial_solution(self):
        """
        Greedy algorithm to get an initial solution (closest-neighbour)        
        """
        cur_node = 0 # start with node 0
        solution = [cur_node]

        free_list = list(self.nodes)
        free_list.remove(cur_node)

        while free_list:
            #print(cur_node)
            closest_dist,cur_node = min([(self.dist_matrix[cur_node][j],j) for j in free_list])
            #closest_dist = min([self.dist_matrix[cur_node][j] for j in free_list])
            #cur_node = self.dist_matrix[cur_node].index(closest_dist)
            free_list.remove(cur_node)
            solution.append(cur_node)

        return solution

    def dist(self,coord1, coord2):
        if self.format_type=='EUC_2D':
            """
            Compute Euclidean distance and round to nearest integer
            """
            dij=int(round(math.sqrt(math.pow(coord1[0] - coord2[0], 2) + math.pow(coord1[1] - coord2[1], 2))))

        elif self.format_type=='GEO':
            """
            Compute GEO Format distance
            """
            # Required constants
            PI = 3.141592
            RRR = 6378.388

            xi=coord1[0]
            xj=coord2[0]
            yi=coord1[1]
            yj=coord2[1]

#            if self.iteration<=1:
#                print("here")
#                print(xi)
#                print(xj)
#                print(yi)
#                print(yj)
#                print()

            # Convert latitude to radians

            deg1 = math.floor(xi)  #multiplying by its sign
            deg2 = math.floor(xj)

            #deg1 = round(xi)
            #deg2 = round(xj)

            mi1 = xi-deg1
            mi2 = xj-deg2

            lati = PI*(deg1+5.0*mi1/3.0)/180.0
            latj = PI*(deg2+5.0*mi2/3.0)/180.0

#            if self.iteration<=1:
#                print(deg1)
#                print(deg2)
#                print(mi1)
#                print(mi2)
#                print(lati)
#                print(latj)
#                print()

            # Convert longitude to radians

            deg1 = round(yi)
            deg2 = round(yj)
            mi1 = yi-deg1
            mi2 = yj-deg2

            #deg1 = (yi/abs(yi))*math.floor(abs(yi))  #multiplying by its sign
            #deg2 = (yj/abs(yj))*math.floor(abs(yj))
            #mi1 = abs(yi-deg1)
            #mi2 = abs(yj-deg2)

            longi = PI*(deg1+5.0*mi1/3.0)/180.0
            longj = PI*(deg2+5.0*mi2/3.0)/180.0

#            if self.iteration<=1:
#                print(deg1)
#                print(deg2)
#                print(mi1)
#                print(mi2)
#                print(longi)
#                print(longj)
#                print()

            # Calculate distance
            q1 = math.cos(longi-longj)
            q2 = math.cos(lati-latj)
            q3 = math.cos(lati+latj)
            dij = int(RRR*math.acos(0.5*((1.0+q1)*q2-(1.0-q1)*q3))+1.0)

#            if self.iteration<=1:
#                print(q1)
#                print(q2)
#                print(q3)
#                print(dij)

        else:
            print('Input data type not currently supported for this project.')
            exit()

        return dij

    def to_dist_matrix(self, coords):
        """
        Returns nxn nested list from a list of length n
        Used as distance matrix: mat[i][j] is the distance between node i and j
        'coords' has the structure [[x1,y1],...[xn,yn]]
        """
        n = len(coords)
        mat = [[self.dist(coords[i], coords[j]) for i in range(n)] for j in range(n)]
        return mat

    def fitness(self, sol):
        """ Objective value of a solution """
        return sum([self.dist_matrix[sol[i - 1]][sol[i]] for i in range(1, self.N)]) + self.dist_matrix[sol[0]][sol[self.N - 1]] #sum of evertyhing in the solution + return back to start node

    def p_accept(self, candidate_fitness):
        """
        Probability of accepting if the candidate is worse than current
        Depends on the current temperature and difference between candidate and current
        """
#        if self.iteration<=10:
#            print(candidate_fitness)
#            print(self.cur_fitness)
#            print(self.T)
#            print(-abs(candidate_fitness - self.cur_fitness) / self.T)
#            print(math.exp(-abs(candidate_fitness - self.cur_fitness) / self.T))
#            print()
        return math.exp(-abs(candidate_fitness - self.cur_fitness) / self.T)

    def accept(self, candidate):
        """
        Accept with probability 1 if candidate is better than current
        Accept with probabilty p_accept(..) if candidate is worse
        """
        candidate_fitness = self.fitness(candidate)
        if candidate_fitness <= self.cur_fitness: # LESS THAN EQUAL
            self.cur_fitness = candidate_fitness
            self.cur_solution = candidate

            if candidate_fitness < self.best_fitness:
                self.best_fitness = candidate_fitness
                self.best_solution = candidate

                # update the time of improves solution
                self.time_fitness_list.append((time.clock()-self.start_time,candidate_fitness))

        else:
            if random.random() < self.p_accept(candidate_fitness):
                self.cur_fitness = candidate_fitness
                self.cur_solution = candidate

    def anneal(self):
        """
        Execute simulated annealing algorithm
        """
        self.start_time=time.clock()
        #for iters in range(self.restarts):
        while (time.clock()-self.start_time)<=self.time_cutoff: # keep running (restarting SA) till time cutoff

            self.T=self.begin_T # start temperature at beginning temperature for every restart
            candidate=list(self.best_solution)#.copy() # initialize every restart with the best solution found so far
            #random.shuffle(candidate)# absolutely random restart
            #print(self.T)

            while self.T >= self.stopping_temperature: # and self.iteration < self.stopping_iter:
                l = random.randint(2, self.N - 1)
                ######## #i = random.randint(0, self.N - l)
                i = random.randint(1, self.N - l)
                candidate[i:(i + l)] = reversed(candidate[i:(i + l)])  # reverse the indexes between l and u
                """
                http://toddwschneider.com/posts/traveling-salesman-with-simulated-annealing-r-and-shiny/

                one way to pick a neighboring tour is to choose two cities on the tour randomly,
                and then reverse the portion of the tour that lies between them.
                This candidate tour might be better or worse compared to the existing tour,
                i.e. shorter or longer.
                """

#                l = random.randint(1, self.N - 1) #CHANGE ORDER OF ALL NODES EXCEPT FIRST
#                i = l
#                while i==l:
#                    i=random.randint(1, self.N - 1) #making sure i is different from l
#
#                l_val=candidate[l]
#                i_val=candidate[i]
#                candidate[l] = i_val #swapping the two nodes to create a neighboring solution
#                candidate[i] = l_val

                self.accept(candidate)
                self.T *= self.alpha
                self.iteration += 1

                self.fitness_list.append(self.cur_fitness)
                candidate = list(self.cur_solution)
            # END OF WHILE LOOP
            self.best_iteration_solution.append(self.best_solution)
            self.best_iteration_fitness.append(self.best_fitness)

        # END OF FOR LOOP / OR While LOOP
        #print('Initial Greedy tour length: ', self.initial_fitness)
        #print('Best tour length: ', self.best_fitness)
        #print('Improvement over greedy solution: ',round((self.initial_fitness - self.best_fitness) / (self.initial_fitness), 4))

def file_reader(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
    coords = []
    for i, line in enumerate(lines):
        #print(i,line)
        if "NAME" in line:
            name = line.split(" ")[-1][:-1]
            continue
        if "DIMENSION" in line:
            dim = int(line.split(" ")[-1])
            continue
        if "EDGE_WEIGHT_TYPE" in line:
            w_type = line.split(" ")[-1][:-1]
        if "NODE_COORD_SECTION" in line:
            for line in lines[i+1:i+1+dim]:
                splits = line.lstrip().split(" ")
                coords.append([float(splits[1]),float(splits[2][:-1])])
            break
    return name,dim,w_type,coords


##################################
############ OUTPUT FILES ########
##################################

def file_output(name,cutoff,rs,best_fitness,best_solution,time_fitness_list):
    ##### Solution File #####
    solFile = '%s_LS1_%s_%s.sol' % (name,cutoff,rs)
    solpath= './OUTPUT/'+solFile
    with open(solpath,'w') as outputFile:
        outputFile.write(str(best_fitness)+'\n')
        outputFile.write(','.join(str(i+1) for i in best_solution)) # adding 1 to adjust for index numbering 0 to 1 etc
    outputFile.close()

    ##### Trace File #####
    traceFile = '%s_LS1_%s_%s.trace' % (name,cutoff,rs)
    tracepath='./OUTPUT/'+traceFile
    with open(tracepath,'w') as outputFile:
        for i in time_fitness_list:
            outputFile.write(str(round(i[0],3))+','+str(i[1])+'\n')
    outputFile.close()

    return None

def LS1(fname, cutoff, rs):
    inp='./DATA/'+fname+".tsp"
    name,dim,w_type,coords=file_reader(inp)

    #t1=time.clock()
    sa = SimAnneal(coords,format_type=w_type,random_seed=rs,time_cutoff=cutoff)
    sa.anneal()
    #print("Total execution time was: ",time.clock()-t1)
    file_output(fname,cutoff,rs,sa.best_fitness,sa.best_solution,sa.time_fitness_list)
    return None