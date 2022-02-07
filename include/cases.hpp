#ifndef FM_TEST_CPP
#define FM_TEST_CPP

#include "fmarch.hpp"

#include <filesystem>
#include <chrono>

namespace fs = std::filesystem;


std::vector<int> getBounds(const std::string& filename, const Vector3D& r0, const Vector3D& v0, const double& dt);
void runAndCompare(const std::string& filename, const std::vector<int> b, const int& it_show, const int& N, const double& vol, const double& b_a);
void meshTest(const int& N, const int& it_show = std::numeric_limits<int>::max());
void flameSpread(const std::string& filename, const pairXY& data, const bool& rem = true);
void tryMesh(const std::string& filemesh, const std::string& fileres, const std::vector<int>& b_in, const functionN& fN, const double& beta = 0);
void tryMesh(const MeshInfo& minfo, const std::string& fileres, const std::vector<int>& b_in, const functionN& fN, const double& beta = 0);
void traverse(const MeshInfo& minfo);
void go2Node(const int& i, const std::vector<Node>& nodes, std::vector<bool>& bb, const uint& depth);
void readFileSolve(const int& argc, char* argv[]);

void performFlameSpreading(const std::string& filename, const std::string& baselog, const std::vector<double>& beta, const std::vector<double>& xs);

#endif