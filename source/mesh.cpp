#include "mesh.hpp"

uint Mesh::isGmsh = 0;

// auxiliary vector for identification of the element that contains a point
const Vector<Vector<int,3>, 4> plane = 
{
    0,1,2,
    0,1,3,
    0,2,3,
    1,2,3
};

/*
Default constructor for the Nodes
*/
Node::Node() : state(0), in_bound(false), coord({0,0,0}), val(std::numeric_limits<double>::infinity()) { }

/*
Constructor the Node object
Takes a vector of strings that contain the nodes information from the
.diff file. Initializes the state at 0 and the value at infinity
Takes the ID, the coordinates and the boundaries from the given string vector
*/
Node::Node(const std::vector<std::string>& f) : Node()
{
    setPars(f);
}

void Node::setPars(const std::vector<std::string>& f)
{
    int n_bound;

    N = std::stoi(f[0])-1;
    // obtains the positional coordinates of the nodes
    for(int i = 2; i < 7; i += 2){
        coord[i/2-1] = std::stod(f[i]);
    }

    
    n_bound = std::stoi(f[8].substr(1, f[8].size()-1));
    in_bound = n_bound;
    for(int i = 0; i < n_bound; ++i){
        bound.push_back(std::stoi(f[9+i]));
    }
}

/*
Connects two nodes (self and other) as being neighbours
*/
void Node::connect(Node& n)
{
    // cannot connect node to self or to an already
    // connected node
    if(N == n.N || findIn(neigh, n.N) != -1) return;

    neigh.push_back(n.N);
    n.neigh.push_back(N);
}

/*
Mesh constructor, takes in a .diff file to read
Will allocate memory for all nodes and elements,
and it will read the file to save the information
about the nodes and meshes, making all necessary
connections 
*/
Mesh::Mesh(const std::string& filename)
{
    std::cout << "Info\t: Reading mesh . . ." << std::endl;
    std::fstream fp(filename, std::ifstream::in);
    std::vector<std::string> aux;
    std::vector<int> aux_n;

    // exits if a bad file is given
    if(fp.fail() || fp.bad()){
        std::cout << "ERROR\t: Could not open file for mesh reading." << std::endl;
        exit(-1);
    }

    // reads the number of elements and nodes
    n_elem = std::stoi(spltString(readLine(fp, 5)).back());
    n_node = std::stoi(spltString(readLine(fp, 0)).back());

    // allocated memory for them
    nodes.resize(n_node, Node());
    elems.reserve(n_elem);

    // skips until the start of the nodes
    while(readLine(fp) != "#");

    for(uint i = 0; i < n_node; ++i){
        // reads, creates and places the node in the vector
        aux = spltString(readLine(fp));
        nodes[std::stoi(aux[0])-1].setPars(aux);
    }

    // skips until the start of the elements
    while(readLine(fp) != "#");

    std::cout << "Info\t: Joining elements . . ." << std::endl;
    for(uint i = 0; i < n_elem; ++i){
        aux = spltString(readLine(fp));
        aux_n.clear();
        for(int j = 3; j < (int)aux.size(); ++j){
            // identifies the nodes in the element
            aux_n.push_back(std::stoi(aux[j])-1);
        }
        // inserts the element
        elems.push_back(aux_n);

        // connects the nodes in the element as being neighbours
        for(int k = 0; k < (int)aux_n.size(); ++k){
            nodes[aux_n[k]].elem.push_back(i);
            for(int l = k+1; l < (int)aux_n.size(); ++l){
                nodes[aux_n[k]].connect(nodes[aux_n[l]]);
            }
        }
    }

    std::cout << "Info\t: Mesh complete." << std::endl; 

    /*
    for(const Node& n : nodes){
        std::cout << "\nNode ID: " << n.N << std::endl;
        std::cout << "In bound: " << n.in_bound << std::endl;
        if(n.in_bound){
            std::cout << "Boundaries: " << n.bound << std::endl;
        }
        std::cout << "In elements: " << n.elem << std::endl;
    }

    std::cout << "\nNumber of elements: " << elems.size() << std::endl; 
    */
}

Mesh::Mesh(const MeshInfo& minfo)
{
    
    std::vector<std::size_t> nTags;
    std::vector<double> coords, pCoords;
    std::vector<int> eTypes, aux, bTags;
    std::vector<std::vector<std::size_t>> eTags, neTags;
    uint pos = 0;
    int btag = -1;
    gmsh::vectorpair dimTags;
    const auto& [f, in_tag, un_tag, mSize] = minfo;
    mSizeFactor = mSize;
    filename = f;

    // initialize Gmsh
    if(isGmsh) gmsh::clear();
    else gmsh::initialize();
    ++isGmsh;
    gmsh::open(f + ".geo");

    // boundary physical groups conditions
    gmsh::model::getEntities(dimTags, 2);
    for(const auto&[dim, tag] : dimTags){
        bTags.emplace_back(tag);
    }
    // if inhibited tag is given, remove surfaces from general boundary
    if(in_tag != -1){
        aux.clear();
        gmsh::model::getEntitiesForPhysicalGroup(2, in_tag, aux);
        for(const int& i : aux){
            removeItem(i, bTags);
        }
    }
    // if initial condition is given, remove surfaces from general boundary
    if(un_tag != -1){
        aux.clear();
        gmsh::model::getEntitiesForPhysicalGroup(2, un_tag, aux);
        for(const int& i : aux){
            removeItem(i, bTags);
        }
    }
    if(bTags.size()){
        btag = gmsh::model::addPhysicalGroup(2, bTags);
    }    

    // meshing operations
    // uses the size callback function to set the global mesh size factor to the given one
    gmm::setSizeCallback([this](int dim, int tag, double x, double y, double z, double lc){ return lc*this->mSizeFactor; });
    // generate 3D mesh
    gmm::generate(3);
    
    // reading nodes
    gmm::getNodes(nTags, coords, pCoords, -1, -1, false, false);
    nodes.resize(nTags.size());
    std::cout << "Info\t: Adding " << nTags.size() << " nodes . . ." << std::endl;
    n_node = nTags.size();
    for(uint i = 0; i < nTags.size(); ++i){
        nodes[nTags[i]-1].N = nTags[i]-1;
        nodes[nTags[i]-1].coord = {coords[3*i], coords[3*i+1], coords[3*i+2]};
    }
    
    // adding boundary conditions to nodes
    if(btag != -1) addBoundary(btag, 0);
    if(in_tag != -1) addBoundary(in_tag, in_tag);
    if(un_tag != -1) addBoundary(un_tag, un_tag);

    // getting elements
    gmm::getElements(eTypes, eTags, neTags);
    for(uint i = 0; i < eTypes.size(); ++i){
        if(eTypes[i] == 4){
            n_elem = eTags[i].size();
            pos = i;
        }
    }

    // joining nodes and adding elements
    elems.reserve(n_elem);
    std::cout << "Info\t: Joining " << n_elem << " elements . . ." << std::endl;
    for(std::size_t i = 0; i < eTags[pos].size(); ++i){
        aux.clear();
        for(std::size_t j = 4*i; j < 4*(i+1); ++j){
            aux.emplace_back(neTags[pos][j]-1);
        }
        elems.emplace_back(aux);
        
        for(int k = 0; k < (int)aux.size(); ++k){
            nodes[aux[k]].elem.push_back(i);
            for(int l = k+1; l < (int)aux.size(); ++l){
                nodes[aux[k]].connect(nodes[aux[l]]);
            }
        }
    }

    std::cout << "Info\t: Done." << std::endl;    
}

/*
Mesh destructor
Closes GMSH if active
*/
Mesh::~Mesh()
{
    if(!--isGmsh) gmsh::finalize();
}


/*
Method that checks in which element a certain point p is located
Can search for both 2D and 3D meshes
If the point is in an edge/face, will return only that edge/face
*/
std::vector<int> Mesh::inElem(const Vector3D& p, const bool& d3D)
{
    bool inElem, inPlane;
    int sc;
    double d;
    Vector3D c, abc;
    std::vector<int> pl_extra;


    if(d3D){
        pl_extra.resize(3);
        for(auto el : elems){
            inElem = true;
            inPlane = false;
            c = {0,0,0};
            for(int e : el){
                c += nodes[e].coord;
                if(nodes[e].coord == p) return {e};
            }
            c = c/el.size();
            for(auto pl : plane){
                abc = cross(nodes[el[pl[1]]].coord - nodes[el[pl[0]]].coord, nodes[el[pl[2]]].coord - nodes[el[pl[0]]].coord);
                d = -(abc*nodes[el[pl[0]]].coord);
                sc = sgn(abc*p + d);
                if(!sc){
                    inPlane = true;
                    for(int i = 0; i < 3; ++i){
                        pl_extra[i] = el[pl[i]];
                    }
                }
                else if(sgn(abc*c + d) != sc){
                    inElem = false;
                    break;
                }
            }

            if(inElem){
                if(inPlane) return pl_extra;
                return el;
            }
        }
    }
    else{
        pl_extra.resize(2);
        for(auto el : elems){
            inElem = true;
            inPlane = false;
            c = {0,0,0};
            for(int e : el){
                c += nodes[e].coord;
                if(nodes[e].coord == p) return {e};
            }
            c = c/el.size();
            for(int i = 0; i < (int)el.size(); ++i){
                abc = nodes[el[i]].coord - nodes[el[(i+1)%el.size()]].coord;
                d = abc[0];
                abc[0] = abc[1];
                abc[1] = -d;
                d = -(abc*nodes[el[i]].coord);
                sc = sgn(abc*p + d);
                if(!sc){
                    inPlane = true;
                    pl_extra = {el[i], el[(i+1)%el.size()]};
                }
                else if(sgn(abc*c + d) != sc){
                    inElem = false;
                    break;
                }
            }

            if(inElem){
                if(inPlane) return pl_extra;
                return el;
            }
        }
    }

    return {};
}

// returns the number of nodes
uint Mesh::getNnode() const
{
	return n_node;
}

// returns the number of elements
uint Mesh::getNelem() const
{
    return n_elem;
}

// resets the mesh
void Mesh::resetMesh()
{
    for(Node& n : nodes){
        n.state = 0;
        n.val = std::numeric_limits<double>::infinity();
    }
}

void Mesh::addBoundary(const int& tag, const int& b)
{
    std::vector<double> coord;
    std::vector<std::size_t> nTags;
    gmm::getNodesForPhysicalGroup(2, tag, nTags, coord);

    for(const std::size_t& n : nTags){
        nodes[n-1].in_bound = true;
        nodes[n-1].bound.emplace_back(b);
    }
}

void Mesh::exportMesh()
{
    if(isGmsh){
        std::cout << "Info\t: Exporting file . . ." << std::endl;
        gmsh::write(filename + ".diff");
        std::cout << "Info\t: Done." << std::endl;
    }
    else std::cout << "Error\t: Mesh was not created. Nothing to export." << std::endl;
}