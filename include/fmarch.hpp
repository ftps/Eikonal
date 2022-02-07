#ifndef FMARCH_HPP
#define FMARCH_HPP

//#include <utility>

//#include <fstream>

#include "mesh.hpp"
#include "heapN.hpp"

// alias of a std::function that receives a coordinate position and a node id, and returns a double
using functionXI = std::function<double(const int&)>;
using functionN = std::function<bool(const Node&)>;
//using tVd = std::tuple<Vector3D, double>;

const functionN fTrue = [](const Node& n){ return true; };

/*
Enumeration of the state each node can have
*/
enum state{far, considered, seen};

/*
Triangulated Fast Marching Method implementation
*/
class FastMarch : public Mesh, private HeapN {
	public:
		FastMarch(const std::string& filename, const functionXI& f, std::vector<int> b, const bool& d3D = true);
		FastMarch(const std::string& filename, const functionXI& f, std::vector<Vector3D> b, const bool& d3D = true);
		FastMarch(const std::string& filename, const std::vector<int>& b_in, const double& beta, const functionN& fN, const bool& d3D = true);
		FastMarch(const MeshInfo& minfo, const std::vector<int>& b_in, const double& beta, const functionN& fN, const bool& d3D = true);
		void reset(const std::vector<int>& b_in, const double& beta, const functionN& fN, const bool& resMesh = true);
		void solve(const int& it_inter = std::numeric_limits<int>::max());
		pairXY getSurface(const int& N);
		double getMaxTime();		
		
	private:
		FastMarch(const std::string& filename, const bool& d3D);
		FastMarch(const MeshInfo& minfo, const bool& d3D);
		void generateHeap(const std::vector<Node*> bound);
		bool updateNode(const int& i, const int& orig = -1);
		double update2P(const std::vector<int>& gnei, const int& i);
		double update3P(const std::vector<int>& gnei, const int& i);
		double length2D(const std::vector<int>& el, const double& t);
		double area3D(const std::vector<int>& el, const double& t);
		//tVd& lookUp(const int& i, const int& j);

		/*
		FMM paramters:
		d3D - 2D or 3D mesh
		f - inverse of the velocity function in the Eikonal equation
		map - look up table for the direction vectors (scrapped)
		max_val - maximum time at the last node, used to compute area curve
		*/
		const bool d3D;
		double beta2;
		functionXI f;
		//std::map<Vector<int,2>, tVd> map;
		double max_val;
		uint bsize;
		std::vector<bool> check;
};

#endif