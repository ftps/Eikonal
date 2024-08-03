#ifndef MESH_HPP
#define MESH_HPP

#include <map>

#include "utils.hpp"
#include "gmsh.h"

//using tVdd = std::tuple<Vector3D, double, double>;
namespace gmm = gmsh::model::mesh;

/*
Data type to represent a node in the mesh for the Eikonal solver
*/
class Node {
    public:
        Node();
        Node(const std::vector<std::string>&);
        void connect(Node&);
        void setPars(const std::vector<std::string>&);

        /*
        Node parameters:
        N - node id
        state - node state as per the enumeration (0 far, 1 considered, 2 seen)
        in_bound - in boundary of the mesh
        neigh - vector containing the id of the neighbouring elements
        elem - vector containing the id of the elements if belongs to
        bound - vector containing the id of the boundaries it belong to (empty if none)
        coord -  3D space coordinates of the node
        val - value of the Eikonal solution attributed to this node
        */
        int N, state;
        bool in_bound;
        std::vector<int> neigh, elem, bound;
        Vector3D coord;
        double val;
};

/*
Struct with information for meshing
file - filename without the ".geo" extension
in_tag - physical group tag of the inhibited surface (-1 and nothing will be inhibited)
un_tag - physical group tag of the initial burning surface (-1 and all non-inhibited surfaces will ignite)
mSize - mesh refinement level as per the global mesh size factor in the UI
*/
struct MeshInfo {
    std::string f = "";
    int in_tag = -1;
    int un_tag = -1;
    double mSize = 0.05;
};

/*
Mesh object that contains all the Nodes and Elements of the mesh
*/
class Mesh {
    public:
        Mesh(const std::string& f);
        Mesh(const MeshInfo& minfo);
        ~Mesh();
        void resetMesh();
        void exportMesh();
        uint getNnode() const;
        uint getNelem() const;
		void exportInfo(std::fstream& fp) const;

        std::vector<int> inElem(const Vector3D& p, const bool& d3D = true);
        //tVdd getVals(const Vector<int,2>& key);

        /*
        Mesh parameters:
        n_elem - number of elements
        n_node - number of nodes
        nodes - vector containing all the node objects
        elems - vector containing vectors that identify the nodes in each element
        map - look up table for the direction vectors and distances between nodes
        isGmsh - indicates if gmsh has been initiased (in case of multiple meshes, future case)
        filename - name in MeshInfo struct
        mSizeFactor - global mesh size factor (as per MeshInfo)
        */
    protected:
        uint n_elem, n_node;
        std::vector<Node> nodes;
        std::vector<std::vector<int>> elems;
    private:
        static uint isGmsh;
        std::string filename;
        double mSizeFactor;

        void addBoundary(const int& tag, const int& b);
};

#endif