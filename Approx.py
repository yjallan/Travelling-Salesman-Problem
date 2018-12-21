'''
CSE6140 Final Project
MST 2-Approximation Implementation
'''

import time
import sys
import math

class RunApproxAlg(object):
    def __init__(self,E,V,cutoff_time,start_time):
        self.edges = E
        self.nodes = V
        self.cutoff = cutoff_time
        self.start_time = start_time
    # Read input data

    # Find minimum spanning tree using Prim's Algorithm
    def computeMST(self):
        # Delete the minimum key from Q
        def removeMin(Q,parent):
            u = min(Q,key = Q.get) # find node with the minimum value
            costu = Q[u]
            par = parent[u]
            del Q[u]
            del parent[u]

            return Q,u,costu,par

        # Initialize dictionary for Q and dictionary for parent
        # Q is updated to reflect current lowest cost neighbor
        # parent is updated to reflect the parent of the lowest cost neighbor
        Q = {}
        parent = {}

        # Add all nodes to Q, with keys = infinity
        for cur in self.nodes:
            Q[cur[0]] = float('inf')

        # Initialize list of already explored nodes
        S = []

        # Set values for first node
        Q[1] = 0
        parent[1] = float('inf') # first node has no parent

        # Initialize MST weight and tree
        MSTweight = 0
        MST = []

        # Run until Q is empty and all nodes are explored
        while len(Q) > 0:
            # Remove the minimum cost node from Q
            [Q,u,costu,par] = removeMin(Q,parent)

            # Add edges (undirectional) to MST
            if len(S) >= 1:
                MST.append((u,par,costu))
                MST.append((par,u,costu))

            # Mark removed node as visited
            S.append(u)

            # Update the cost of the MST
            MSTweight = MSTweight+costu

            # Find neighbors of node that was just explored
            for v in self.edges[u]:
                neighbor = v[0]
                if neighbor not in S:
                    # Get the cost of the edge from u to v
                    cost = v[1] # v is (node, cost)

                    # Update the value of the neighbor in Q if it is less than the existing value
                    if cost < Q[neighbor]:
                        Q[neighbor] = cost
                        parent[neighbor] = u

        return MST

    # Christofides algorithm
    def christofides(self,T):
        # Get the cost of an edge
        def GetCost(u,v):
            for i in self.edges[u]:
                if i[0] == v:
                    edgeCost = i[1]
            return edgeCost

        # Traverse a graph using depth first search given a start node
        # Construct the final tour and calculate its length
        def DFS(start):
            S = []
            S.append(start)
            visited = []
            foundTour = []
            foundCost = 0
            while S:
                cur = S.pop()
                if cur not in visited:
                    visited.append(cur)
                    if len(visited) == 1:
                        prev = cur
                    else: # If this is not the first node added
                        curCost = GetCost(prev,cur)
                        foundTour.append((prev,cur,curCost))
                        foundCost += curCost
                        prev = cur

                    sortedges = dictT[cur]
                    sortedges.sort(key = lambda a: a[1])
                    for i in reversed(sortedges):
                        S.append(i[0])

            # Close the cycle
            lastCost = GetCost(prev,start)
            foundCost = foundCost + lastCost
            foundTour.append((prev,start,lastCost))
            return visited,foundCost

        # dictionary of all edges in the MST (each edge counted twice, once for each direction)
        dictT = dict((i+1,[]) for i in range(len(self.edges)))
        for edge in T:
            if edge[0] in dictT:
                dictT[edge[0]].append((edge[1],edge[2]))
            elif edge[1] in dictT:
                dictT[edge[1]].append((edge[0],edge[2]))

        # Find nodes with odd degree
        oddDeg = []
        for node in dictT:
            curDeg = len(dictT[node])
            if (curDeg % 2) != 0:
                oddDeg.append(node)

        # Get edges without repeats for odd tree of MST
        oddEdges = []
        for node in oddDeg:
            for v in self.edges[node]:
                if v[0] in oddDeg and (v[0], node, v[1]) not in oddEdges:
                    oddEdges.append((node,v[0],v[1]))

        # Sort edges in increasing order by edge weight for finding minimum matching
        oddEdges.sort(key = lambda a: a[2])

        # Set an initial minimum matching value at infinity
        mintest = float('inf')

        # Find a minimum matching of a subgraph of the initial with only the odd degree nodes in the MST
        for each in oddEdges:
            # Initialize variables
            found = []
            found.append(each[0])
            found.append(each[1])
            count = 1
            matching = []
            matching.append(each)
            newtest = each[2]

            # Find sums of 1/2*number-of-odd-degree-nodes weights to find minimum matching
            for test in oddEdges:
                # If the nodes have not yet been added to the matching and the perfect match size has not been reached
                if test[0] not in found and test[1] not in found and count <= len(oddDeg)/2:
                    count += 1
                    found.append(test[0])
                    found.append(test[1])
                    newtest = newtest + test[2]
                    matching.append(test)

                    # If no more edges can be added to the matching, update the new minimum if needed
                    if newtest < mintest and len(matching) == len(oddDeg)/2:
                        mintest = newtest
                        minmatch = matching
                        # If the time limit has run out, continue on to calculating the final tour
                        if (time.time()-self.start_time) > self.cutoff:
                            break

        # Add the minimum match edges to the MST
        for edge in minmatch:
            if edge[0] in dictT:
                addedge = GetCost(edge[0],edge[1])
                dictT[edge[0]].append((edge[1],addedge))

        # Initialize final tour and tour length
        tourcost = float('inf') # initial cost is 0
        tested = []

        # Create tours starting at each node
        # Run DFS from each to determine which produces shortest cycle
        for root in dictT:
            if root not in tested:
                tested.append(root)
                tour,fincost = DFS(root)

                # Update the tour if starting from this root was better than the others
                if fincost < tourcost:
                    tourcost = fincost
                    visited = tour

            # If the time limit has run out, return the last tour calculated
            if (time.time()-self.start_time) > self.cutoff:
                break

        # return final tour and weight
        return visited, tourcost

def readData(filename):
    filename = "DATA/" + filename + ".tsp"
    with open(filename,'r') as readfile:
        # Initialize variables for nodes and edges
        # nodes: list, edges: dictionary
        nodes = []
        edges = {}
        for line in readfile:
            line = line.strip('\n').split(' ')
            if line[0] == '':
                line.pop(0)
            # Indicate type of distance calculation
            if line:
                if line[0] == 'EDGE_WEIGHT_TYPE:':
                    ewt = line[1]
                # Save nodes and add a key for each node in edges
                elif line[0].isdigit():
                    cur = int(line[0])
                    x = float(line[1])
                    y = float(line[2])
                    nodes.append((cur,x,y))
                    edges[cur] = []

    # Calculate distances
    def calc_dist(xi,xj,yi,yj,weight_type):
        # Compute Euclidean distance and round to nearest integer
        if weight_type == 'EUC_2D':
            xd = xi-xj
            yd = yi-yj
            dij = int(round(math.sqrt(xd*xd +yd*yd)))
        # Compute geographical distance for coordinates given in lat/long coordinate pairs
        elif weight_type == 'GEO':
            # Required constants
            PI = 3.141592
            RRR = 6378.388

            # Convert latitude to radians
            deg1 = math.floor(xi)
            deg2 = math.floor(xj)
            mi1 = xi-deg1
            mi2 = xj-deg2
            lati = PI*(deg1+5.0*mi1/3.0)/180.0
            latj = PI*(deg2+5.0*mi2/3.0)/180.0

            # Convert longitude to radians
            deg1 = round(yi)
            deg2 = round(yj)
            mi1 = yi-deg1
            mi2 = yj-deg2
            longi = PI*(deg1+5.0*mi1/3.0)/180.0
            longj = PI*(deg2+5.0*mi2/3.0)/180.0

            # Calculate distance
            q1 = math.cos(longi-longj)
            q2 = math.cos(lati-latj)
            q3 = math.cos(lati+latj)
            dij = int(RRR*math.acos(0.5*((1.0+q1)*q2-(1.0-q1)*q3))+1.0)

        # If data is provided that is not in EUC_2D or GEO format...
        else:
            print('Input data type not currently supported for this project.')
            exit()

        return dij

    # Add all edges (both directions) to dictionary
    # u : (v,distance)
    for i in nodes:
        for j in nodes:
            cost = calc_dist(i[1],j[1],i[2],j[2],ewt)
            if i!=j:
                edges[i[0]].append((j[0],cost))

    return edges, nodes

def writeOutput(tour,tcost,timeResults,instance,cutoff):
    # Write results to output file
    solFile = 'OUTPUT/'+'%s_Approx_%s.sol' % (instance,cutoff)
    with open(solFile,'w') as outputFile:
        outputFile.write(str(tcost)+'\n')
        outputFile.write(','.join(str(i) for i in tour))

    traceFile = 'OUTPUT/'+'%s_Approx_%s.trace' % (instance,cutoff)
    with open(traceFile,'w') as outputFile:
        outputFile.write('{:1.15f}'.format(timeResults)+'\n')
        outputFile.write('{:1.15f}'.format(timeResults)+','+str(tcost))


def Approx(fname, cutoff):
    # Find approximate TSP solution
    edges,nodes = readData(fname)

    # Time the algorithm
    startTime = time.time()
    approx = RunApproxAlg(edges,nodes,cutoff,startTime)
    MST = approx.computeMST()
    tour,tcost = approx.christofides(MST)
    totalTime = (time.time()-startTime)

    # Write outputs
    instanceName = fname.split('.')[0]
    writeOutput(tour,tcost,totalTime,instanceName,cutoff)

# Approx('NYC.tsp',600)
