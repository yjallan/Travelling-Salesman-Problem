from sys import argv
import time
import math
import random
import os

class Node:
	"""
	represents a node in a TSP tour
	"""
	def __init__(self, coords):
		self.num = coords[0]-1 # start position in a route's order
		self.x = coords[1]   # x coordinate
		self.y = coords[2]   # y coordinate

	def __str__(self):
		"""
		returns the string representation of a Node
		"""
		return str(self.num)

	def __eq__(self, other):
		return self.__dict__ == other.__dict__

	def euclidean_dist(self, other):
		"""
		returns the Euclidean distance between this Node and other Node
		other - other node
		"""
		dx = self.x - other.x
		dy = self.y - other.y
		return int(round((dx**2+dy**2)**.5))

	def geo_dist(self, other):

		# Required constants
	    PI = 3.141592
	    RRR = 6378.388
	    xi, xj, yi, yj=self.x, other.x, self.y, other.y

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
	    return int(RRR*math.acos(0.5*((1.0+q1)*q2-(1.0-q1)*q3))+1.0)

def parse_input_route(filename, SEED):
	"""
	returns initial route as read from input file, None if parsing errors occur
	filename - name of the input file with '.tsp' extension
	"""
	f = open(filename, 'r')
	route = []
	dimension = -1 
	dimension_found = False
	node_section_found = False

	# Parse header
	for line in f:
		if "DIMENSION" in line:
			tokens = line.split()
			dimension = int(tokens[-1])
			dimension_found = True
		if "EDGE_WEIGHT_TYPE" in line:
			tokens = line.split()
			weight = tokens[-1]
		if "NODE_COORD_SECTION" in line:
			node_section_found = True			
			break

	# Check for parsing errors in header
	if not dimension_found:
		print("Parsing error: DIMENSION not found")
		f.close()
		return None
	elif not node_section_found:
		print("Parsing error: NODE_COORD_SECTION header not found")
		f.close()
		return None

	# Parse nodes
	for line in f:
		if "EOF" in line:
			break
		coords = get_coords(line)
		if not coords:
			print("Parsing error: Invalid node data found")
			f.close()
			return None
		route.append(Node(coords))
	f.close()

	# Check for parsing error with nodes
	if len(route) != dimension:
		print("Parsing error: number of nodes found does not match dimension")
		return None
	random.seed(SEED)
	random.shuffle(route)
	return route, weight

def get_coords(line):
	"""
	returns the line data as numerals, None if line contains more than 
		3 items or non-numerics in the line
	line - string containing the data
	"""
	data = line.split()
	if len(data) == 3:
		try:
			coords = (int(data[0]), int(float(data[1])), int(float(data[2])))
			return coords
		except ValueError:
			pass
	return None

def route_distance(route, weight):
	"""
	returns the distance traveled for a given tour
	route - sequence of nodes traveled, does not include
	        start node at the end of the route
	"""
	dist = 0
	prev = route[-1]
	if weight=='EUC_2D':
		for node in route:
			dist += node.euclidean_dist(prev)
			prev = node
		return dist
	else:
		for node in route:
			dist += node.geo_dist(prev)
			prev = node
		return dist		

def swap_2opt(route, i, k):
	"""
	swaps the endpoints of two edges by reversing a section of nodes, 
		ideally to eliminate crossovers
	returns the new route created with a the 2-opt swap
	route - route to apply 2-opt
	i - start index of the portion of the route to be reversed
	k - index of last node in portion of route to be reversed
	pre: 0 <= i < (len(route) - 1) and i < k < len(route)
	post: length of the new route must match length of the given route 
	"""
	assert i >= 0 and i < (len(route) - 1)
	assert k > i and k < len(route)
	new_route = route[0:i]
	new_route.extend(reversed(route[i:k + 1]))
	new_route.extend(route[k+1:])
	assert len(new_route) == len(route)
	return new_route

def run_2opt(route, weight, tracefile, start, cutoff):
	"""
	improves an existing route using the 2-opt swap until no improved route is found
	best path found will differ depending of the start node of the list of nodes
		representing the input tour
	returns the best path found
	route - route to improve
	"""
	with open(tracefile, 'w') as file:
		improvement = True
		best_route = route
		best_distance = route_distance(route, weight)
		while improvement and time.clock()-start<=cutoff: 
			improvement = False
			for i in range(len(best_route) - 1):
				for k in range(i+1, len(best_route)):
					new_route = swap_2opt(best_route, i, k)
					new_distance = route_distance(new_route, weight)
					if new_distance < best_distance:
						best_distance = new_distance
						best_route = new_route
						file.write(str(round(time.clock()-start,6))+', '+str(new_distance)+'\n')
						improvement = True
						break #improvement found, return to the top of the while loop
				if improvement:
					break
		assert len(best_route) == len(route)
	return best_route

def print_results(route, outputfile, weight):
	"""
	write into output file the nodes in the final route and route information
	route - route to print
	outputfile - name of the original input filename
	weight - EUC_2D or GEO

	"""
	with open(outputfile, 'w') as file:
		file.write(str(route_distance(route, weight))+'\n')
		for node in route[:-1]:
			file.write(str(node.num+1)+', ')
		file.write(str(route[-1].num+1))
		

def LS2(fname, cutoff, seed):
	inp='./DATA/'+fname+'.tsp'
	tracefile = './OUTPUT/'+fname+'_LS2_'+str(cutoff)+'_'+str(seed)+'.trace'
	solfile = './OUTPUT/'+fname+'_LS2_'+str(cutoff)+'_'+str(seed)+'.sol'
	route, weight = parse_input_route(inp,seed)

	# Run 2opt
	start = time.clock() #start time of running 2opt
	route = run_2opt(route, weight, tracefile, start, cutoff)
	end = time.clock()   #end time of running 2opt
	print_results(route, solfile, weight)
