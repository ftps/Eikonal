#include "fmarch.hpp"

/*
Base private constructors that creates the mesh and the heap
*/
FastMarch::FastMarch(const std::string& filename, const bool& d3D) : Mesh(filename), HeapN([this](const int& i){ return nodes[i].val; }), d3D(d3D) {}
FastMarch::FastMarch(const MeshInfo& minfo, const bool& d3D) : Mesh(minfo), HeapN([this](const int& i){ return nodes[i].val; }), d3D(d3D) {}

/*
First FMM constructor that takes in the filename of the mesh,
the inverse of the velocity and a list of the boundaries that
are initially ignited
The constructor finds all points in the given boundaries, sets
them as seen  and sends them to the heap generator
*/
FastMarch::FastMarch(const std::string& filename, const functionXI& f, std::vector<int> b, const bool& d3D) : FastMarch(filename, d3D)
{
	std::vector<Node*> bound;

	this->f = f;
	std::cout << "Info\t: Creating initial boundary . . ." << std::endl;
	for(Node& n : nodes){
		if(n.in_bound){
			for(int bb : n.bound){
				if(findIn(b, bb) != -1){
					n.val = 0;
					n.state = seen;
					bound.push_back(&n);
				}
			}
		}
	}

	generateHeap(bound);
}

/*
Second FMM constructor that takes in the filename of the mesh,
the inverse of the velocity and a list of point locations that
are initially ignited
The constructor finds the elements that contain the given points,
sets them as seen and interpolates their value, and sends them
to the heap generator
*/
FastMarch::FastMarch(const std::string& filename, const functionXI& f, std::vector<Vector3D> b, const bool& d3D) : FastMarch(filename, d3D)
{
	std::vector<Node*> bound;

	this->f = f;
	std::cout << "Info\t: Creating initial boundary . . ." << std::endl;
	for(Vector3D& v : b){
		for(int i : inElem(v, d3D)){
			nodes[i].state = seen;
			nodes[i].val = f(i)*modV(nodes[i].coord - v);
			bound.push_back(&(nodes[i]));
		}
	}

	generateHeap(bound);
}

FastMarch::FastMarch(const std::string& filename, const std::vector<int>& b_in, const double& beta, const functionN& fN, const bool& d3D) : FastMarch(filename, d3D)
{
	reset(b_in, beta, fN, false);
}

FastMarch::FastMarch(const MeshInfo& minfo, const std::vector<int>& b_in, const double& beta, const functionN& fN, const bool& d3D) : FastMarch(minfo, d3D)
{
	reset(b_in, beta, fN, false);
}

/*
Heap generator method
This method creates the heap, takes the list of all boundary points and
starts interpolating the values of their neighbours, adding them to the heap
*/
void FastMarch::generateHeap(const std::vector<Node*> bound)
{
	functionW w = [this](int i){return this->nodes[i].val;};
	std::cout << "Info\t: Initialising heap . . ." << std::endl;
	reserve(nodes.size());

	std::cout << "Info\t: Adding boundary to heap and updating neighbours . . ." << std::endl;
	std::cout << "Info\t: Boundary size: " << bound.size() << std::endl;
	std::cout << "Info\t: Total mesh size: " << nodes.size() << std::endl;
	bsize = bound.size();
	for(const Node* n : bound){
		check[n->N] = true;
		for(int i : n->neigh){
			if(nodes[i].state == far){
				if(updateNode(i)){
					addNew(i);
				}
			}
		}
	}
}

/*
Main solver for the FMM
The solver extracts the root of the priority key as the best approximation
and then updates the surrounding nodes and if successful, updates their
position in the heap or adds them as a new value
Prints the current iteration at it_inter intervals, if no value is specified
it will never print the iterations
*/
void FastMarch::solve(const int& it_inter)
{
	int a, i = 0;

	std::cout << "Info\t: Starting solver . . ." << std::endl;
	while(!isEmpty()){
		++i;
		a = extractRoot(); // extracts best possible approximation
		nodes[a].state = seen; // sets its state as seen
		check[a] = true;
		if(!(i%it_inter)){
			std::cout << "Info\t: In iteration " << i << " of a max " << n_node-bsize << std::endl;
		}
		for(int k : nodes[a].neigh){	// updates its neighbours
			if(nodes[k].state == seen) continue;
			else if(updateNode(k, a)) addNew(k);
			else if(nodes[k].state == considered) siftUp(k);
		}
	}
	std::cout << "Info\t: Final iteration: " << i << std::endl << std::endl;
	max_val = nodes[a].val;

	a = 0;
	i = 0;
	for(uint j = 0; j < n_node; ++j){
		if(!check[j]){
			++a;
			if(nodes[j].elem.size() == 0) ++i;
		}
	}

	if(a){
		std::cout << "Warning\t: A total of " << a << " nodes have not been visited." << std::endl;
		std::cout << "Warning\t: " << i << " were not part of an element." << std::endl;
	}
	
}

/*
Method that updates a node i
If the node orig is defined, it will only consider updates
that contain the orig node in the calculations (used to
reduce computational necessity of the algoritm during
updates). 
*/
bool FastMarch::updateNode(const int& i, const int& orig)
{
	bool neue = false, in_orig, s;
	double aux;
	std::vector<double> upds;
	std::vector<int> gnei;

	// gives information if the nodes has already been inserted into the queue or not
	if(nodes[i].state == far){
		nodes[i].state = considered;
		neue = true;
	}

	// updates according to each element the node belongs to
	for(int ele : nodes[i].elem){
		gnei.clear();
		in_orig = false;
		for(int j : elems[ele]){
			if(nodes[j].state == seen){
				if(orig == j) in_orig = true;
				gnei.push_back(j);
			}
		}
		// if no other nodes in the element are seen or the element doesnt include the newly seen node, ignores
		if(!(orig == -1 || in_orig) || !gnei.size()) continue;
		else if(gnei.size() == 1){ // one point update
			upds.push_back(nodes[gnei[0]].val + f(i)*modV(nodes[i].coord - nodes[gnei[0]].coord));
			continue;
		}
		else if(gnei.size() == 2) aux = update2P(gnei, i); // two point update
		else aux = update3P(gnei, i); // three point update

		if(aux >= 0) upds.push_back(aux); // if the update is successful, adds it to the possible updates

	}

	// equal to 's, aux = min(upds)' in python, struct unpacking
	std::tie(s, aux) = min(upds);
	if(!s){
		// checks if there is no update and if the node is new, cancels the insertion
		if(neue){
			nodes[i].state = far;
			neue = false;
		}
	}
	else{
		// changes the value of the node to the smallest of the updates if it's smaller than the current value
		nodes[i].val = (aux < nodes[i].val) ? aux : nodes[i].val;
	}

	return neue;
}

/*
Two point update of a node
Returns -1 if a valid update isn't possible
*/
double FastMarch::update2P(const std::vector<int>& gnei, const int& i)
{
	Vector3D P;
	Vector<Vector3D, 2> PP; 
	Vector<double, 2> a, b;
	Matrix<double, 2> Q;
	double h, c1, c2, c3;

	// generates the vectors a and b, and the unitary direction vectors

	P = nodes[i].coord - nodes[gnei[0]].coord;
	h = modV(P);
	a[0] = 1/h;
	b[0] = -nodes[gnei[0]].val/h;
	PP[0] = P/h;

	P = nodes[i].coord - nodes[gnei[1]].coord;
	h = modV(P);
	a[1] = 1/h;
	b[1] = -nodes[gnei[1]].val/h;
	PP[1] = P/h;

	// generates matrix Q^-1
	Q[0][0] = 1;
	Q[0][1] = PP[0]*PP[1];
	Q[1][0] = Q[0][1];
	Q[1][1] = 1;

	// inverts to get Q
	Q = inverse(Q);

	// Eikonal quadratic coefficients
	c1 = a*(Q*a);
	c2 = a*(Q*b);
	c3 = b*(Q*b) - sqr(f(nodes[i].N));

	// quadratic determinant, checks for real solution
	h = sqr(c2) - c1*c3;
	if(h < 0) return -1;
	h = (-c2 + sqrt(h))/c1;
	
	// upwind criterion
	for(double p : Q*(a*h + b)){
		if(p < 0) return -1;
	}

	return h;
}

/*
Three point update of a node
Returns -1 if a valid update isn't possible
*/
double FastMarch::update3P(const std::vector<int>& gnei, const int& i)
{
	Vector<Vector3D, 3> PP; 
	Vector3D a, b, P;
	Matrix<double, 3> Q;
	double h, c1, c2, c3;

	// generates the vectors a and b, and the unitary direction vectors
	for(int k = 0; k < 3; ++k){
		P = nodes[i].coord - nodes[gnei[k]].coord;
		h = modV(P);
		a[k] = 1/h;
		b[k] = -nodes[gnei[k]].val/h;
		PP[k] = P/h;
	}

	// generates the matrix Q^-1
	for(int k = 0; k < 3; ++k){
		for(int l = 0; l < 3; ++l){
			Q[k][l] = PP[k]*PP[l];
		}
	}

	// inverts to get Q
	Q = inverse(Q);

	// Eikonal quadratic coefficients
	c1 = a*(Q*a);
	c2 = a*(Q*b);
	c3 = b*(Q*b) - sqr(f(nodes[i].N));

	// quadratic determinant, checks for real solution
	h = sqr(c2) - c1*c3;
	if(h < 0) return -1;
	h = (-c2 + sqrt(h))/c1;

	// upwind criterion
	for(double p : Q*(a*h + b)){
		if(p < 0) return -1;
	}

	return h;
}


/*
Function a map of the surface area as a function of time at
N equally spaced points. It does this by moving throught every
element in the mesh and calculating the cross-section of the
curve at each point in time that the surface intersects that
element. It does this for both 2D (length) and 3D elements
*/
pairXY FastMarch::getSurface(const int& N)
{
	double dt = max_val/(N-1);
	std::vector<double> minmax;
	pairXY res;
	int i_min, i_max;

	// initializes the points with the correct time point and area as 0
	for(int i = 0; i < N; ++i){
		res.emplace_back(i*dt, 0);
	}

	std::cout << "Info\t: Generating burn area data . . ." << std::endl;

	for(const std::vector<int>& el : elems){
		minmax.clear();
		// gets time values at each node
		for(int n : el){
			minmax.push_back(nodes[n].val);
		}

		// gets minimum and maximum time value at this element
		auto [b1, min_d] = min(minmax);
		if(!b1) exit(-2);
		auto [b2, max_d] = max(minmax);
		if(!b2) exit(-2);

		// gets minimum and maximum time points
		i_min = (int)(floor((N-2)*min_d/max_val));
		i_max = (int)(floor((N-2)*max_d/max_val))+2;

		// increments the area at each time point with the
		// corresponding value at that element
		for(int i = i_min; i < i_max; ++i){
			if(d3D){
				res[i].second += area3D(el, res[i].first);
			}
			else{
				res[i].second += length2D(el, res[i].first);
			}
		}
	}

	std::cout << "Info\t: Done." << std::endl;
	return res;
}

/*
Calculates the length of the curve intersecting an element at time t (2D case). 
*/
double FastMarch::length2D(const std::vector<int>& el, const double& t)
{
	std::vector<int> up, bl, *two;
	Vector3D d1, d2;
	int one;

	for(int n : el){
		if(nodes[n].val > t) up.push_back(n);
		else bl.push_back(n);
	}

	if(!(bl.size() && up.size())) return 0;
	if(up.size() == 2){
		two = &up;
		one = bl[0];
	}
	else{
		two = &bl;
		one = up[0];
	}

	d1 = (nodes[(*two)[0]].coord - nodes[one].coord)*((t - nodes[one].val)/(nodes[(*two)[0]].val - nodes[one].val));
	d2 = (nodes[(*two)[1]].coord - nodes[one].coord)*((t - nodes[one].val)/(nodes[(*two)[1]].val - nodes[one].val));

	return modV(d2 - d1);
}

/*
Calculates the area of the surface intersecting an element that at time t.
*/
double FastMarch::area3D(const std::vector<int>& el, const double& t)
{
	std::vector<int> up, bl, *three;
	int one;
	Vector<Vector3D, 3> p;
	Matrix<Vector3D, 2> P;

	// spliting the nodes into in front of behind the wave-front
	for(int n : el){
		if(nodes[n].val > t) up.push_back(n);
		else bl.push_back(n);
	}

	// if there are no points in front or behind the wave-front, no intersection
	if(!(bl.size() && up.size())) return 0;
	// two points on each side, 4 edge case
	else if(bl.size() == up.size()){
		for(int i = 0; i < 2; ++i){
			for(int j = 0; j < 2; ++j){
				// inerpolation of wave-front location in the edges
				P[i][j] = (nodes[up[i]].coord - nodes[bl[j]].coord)*((t - nodes[bl[j]].val)/(nodes[up[i]].val  - nodes[bl[j]].val)) + nodes[bl[j]].coord;			}
		}

		// area calculation using cross product
		return 0.5*(modV(cross(P[0][1] - P[0][0], P[1][0] - P[0][0])) + modV(cross(P[0][1] - P[1][1], P[1][0] - P[1][1])));
	}

	// remaining case, 3 edge case
	// selection of single and triple points location
	if(up.size() == 3){
		three = &up;
		one = bl[0];
	}
	else{
		three = &bl;
		one = up[0];
	}

	for(int i = 0; i < 3; ++i){
		// inerpolation of wave-front location in the edges
		p[i] = (nodes[(*three)[i]].coord - nodes[one].coord)*((t - nodes[one].val)/(nodes[(*three)[i]].val - nodes[one].val));
	}

	// area calculation using cross product
	return 0.5*modV(cross(p[2] - p[0], p[1] - p[0]));
}

double FastMarch::getMaxTime()
{
	return max_val;
}

void FastMarch::reset(const std::vector<int>& b_in, const double& beta, const functionN& fN, const bool& resMesh)
{
	if(resMesh){
		std::cout << "Info\t: Clearing mesh . . ." << std::endl;
		resetMesh();
		for(uint i = 0; i < check.size(); ++i) check[i] = false;
		while(!isEmpty()) extractRoot();
		
	}

	std::vector<Node*> bound;
	beta2 = 1.0/(1.0 + beta);
	check.resize(n_node, false);

	f = [this, &b_in](const int& i){
		bool not_in = !this->nodes[i].in_bound;
		if(this->nodes[i].in_bound){
			for(int b : this->nodes[i].bound){
				if(findIn(b_in, b) != -1){not_in = true; break;}
			}
		}
		return (not_in) ? 1.0 : this->beta2;
	};

	std::cout << "Info\t: Creating initial boundary . . ." << std::endl;
	for(Node& n : nodes){
		if(n.in_bound){
			for(int bb : n.bound){
				if(findIn(b_in, bb) == -1 && fN(n)){
					n.val = 0;
					n.state = seen;
					bound.push_back(&n);
					break;
				}
			}
		}
	}

	generateHeap(bound);
}








/*
Function created as a look up table to obtain the direction vectors
and node distances to void repeating calculation. Didn't result in
a performance increase so the idea was scrapped.
*//*
tVd& FastMarch::lookUp(const int& i, const int& j)
{
	Vector<int,2> key = {i,j};
    auto val = map.find(key);
    Vector3D V;
    double h1;

    if(val != map.end()){
        return val->second;
    }

    V = nodes[j].coord - nodes[i].coord;
    h1 = 1/modV(V);

    map[key] = std::make_tuple(V*h1, h1);

    return map[key];
}
*/

