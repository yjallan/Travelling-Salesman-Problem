'''
CSE6140 Final Project
Exact Solution - Branch and Bound Implementation
'''

import numpy as np
import math
import time

def file_reader(fname):


    fname = "DATA/" + fname + ".tsp"

    with open(fname, 'r') as f:
        lines = f.readlines()
    x = []
    y = []
    for i, line in enumerate(lines):
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
                splits = line.split(" ")
                x.append(float(splits[1]))
                y.append(float(splits[2][:-1]))
            break
    return name, dim, w_type, x, y


def EUC_2D(x1, y1, x2, y2):
    xd = x1 - x2
    yd = y1 - y2
    d12 = round(math.sqrt(xd**2 + yd**2))
    return int(d12)


def GEO(xi, yi, xj, yj):
    PI = 3.141592
    RRR = 6378.388

    # Convert latitude to radians
    deg1 = math.floor(xi)
    deg2 = math.floor(xj)
    mi1 = xi - deg1
    mi2 = xj - deg2
    lati = PI * (deg1 + 5.0 * mi1 / 3.0) / 180.0
    latj = PI * (deg2 + 5.0 * mi2 / 3.0) / 180.0

    # Convert longitude to radians
    deg1 = round(yi)
    deg2 = round(yj)
    mi1 = yi-deg1
    mi2 = yj-deg2
    longi = PI * (deg1 + 5.0 * mi1 / 3.0) / 180.0
    longj = PI * (deg2 + 5.0 * mi2 / 3.0) / 180.0

    # Calculate distance
    q1 = math.cos(longi-longj)
    q2 = math.cos(lati-latj)
    q3 = math.cos(lati+latj)
    dij = int(RRR * math.acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0)

    return dij


def Adj(coords, w_type):
    if w_type == "EUC_2D":
        func = EUC_2D
    else:
        func = GEO
    dim  = len(coords)
    adj = np.zeros((dim, dim)).astype(np.int)

    for i in range(dim-1):

        adj[i, i] = 1e8

        x1, y1 = coords[i]

        for j in range(i+1, dim):
            x2, y2 = coords[j]
            adj[i, j] = func(x1, y1, x2, y2)
            adj[j, i] = adj[i, j]

    adj[dim-1, dim-1] = 1e8

    return adj


def firstMin(adj, i, N, imax = 1e16):
    minval = imax
    for k in range(N):
        if (adj[i,k] < minval) and (i != k):
            minval = adj[i,k];
    return minval


def secondMin(adj, i, N, imax = 1e16):
    first = imax
    second = imax
    for j in range(N):
        if (i == j):
            continue
        if (adj[i,j] <= first):
            second = first
            first = adj[i,j]
        elif (adj[i,j] <= second) and (adj[i,j] != first):
            second = adj[i,j]
    return second

def writeOutput(tour, tcost, timeResults, instance,cutoff):
    # Write results to output file
    tour = [i+1 for i in tour]

    solFile = '%s_BnB_%s.sol' % (instance,cutoff)
    with open('./OUTPUT/'+solFile,'w') as outputFile:
        outputFile.write(str(tcost)+'\n')
        outputFile.write(','.join(str(i) for i in tour))

    traceFile = '%s_BnB_%s.trace' % (instance,cutoff)
    with open('./OUTPUT/'+traceFile,'w') as outputFile:
        outputFile.write('{:1.15f}'.format(timeResults)+'\n')
        outputFile.write('{:1.15f}'.format(timeResults)+','+str(tcost))

def BranchAndBound(fname, cutoff=10):
    name, dim, w_type, x, y = file_reader(fname)
    x, y = np.asarray(x)[:,None], np.asarray(y)[:,None]
    coords = np.concatenate([x, y], axis=1)
    adj = Adj(coords, w_type)
    N = adj.shape[0]

    bnb = BnB(adj, N, cutoff)

    bnb.TSP()
    tour = bnb.tour
    cost = bnb.cost
    writeOutput(tour, cost, bnb.elpsd, fname, cutoff)



class BnB():
    def __init__(self, adj, N, cutoff=10):
        self.adj = adj
        self.tour = []
        self.N = N
        self.cost = float('inf')
        self.visited = set()
        self.cutoff_time = cutoff
        self.elpsd = 0


    def TSP_rec(self, curr_bound, curr_weight, level, curr_path, strt_time):
        if level == self.N or (time.clock()-strt_time) >= self.cutoff_time:
            curr_cost = curr_weight + self.adj[curr_path[self.N-1], curr_path[0]]
            if curr_cost < self.cost:
                self.cost = curr_cost
                self.tour = curr_path.copy()
            return

        for i in range(self.N):
            if i not in self.visited:
                temp = curr_bound
                curr_weight += self.adj[curr_path[level-1]][i]
                if level == 1:
                    curr_bound -= (firstMin(self.adj, curr_path[0], self.N) + firstMin(self.adj, i, self.N)) * 0.5
                else:
                    curr_bound -= (secondMin(self.adj, curr_path[level-1], self.N) + firstMin(self.adj, i, self.N)) * 0.5

                if (curr_bound + curr_weight) < self.cost:
                    curr_path[level] = i
                    self.visited.add(i)
                    self.TSP_rec(curr_bound, curr_weight, level+1, curr_path, strt_time)

                curr_weight -= self.adj[curr_path[level-1], i]
                curr_bound = temp
                self.visited = set(curr_path[:level])


    def TSP(self):

        strt_time = time.clock()
        curr_path = [-1] * (self.N + 1)
        curr_bound = 0

        for i in range(self.N):
            curr_bound += firstMin(self.adj, i, self.N) + secondMin(self.adj, i, self.N)

        curr_bound *= 0.5
        curr_path[0] = 0
        self.visited.add(0)
        self.TSP_rec(curr_bound, 0, 1, curr_path, strt_time)
        self.tour[-1] = (0)
        self.elpsd = time.clock() - strt_time
