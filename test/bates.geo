// Gmsh project created on Mon Apr  8 23:23:51 2024
SetFactory("OpenCASCADE");
//+
// Dimensions of BATES grain
R = 2;
r = 0.5;
L = 6;
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, r, 0, 1.0};
//+
Point(3) = {0, -r, 0, 1.0};
//+
Point(4) = {0, 0, r, 1.0};
//+
Point(5) = {0, 0, -r, 1.0};
//+
Point(6) = {0, R, 0, 1.0};
//+
Point(7) = {0, -R, 0, 1.0};
//+
Point(8) = {0, 0, R, 1.0};
//+
Point(9) = {0, 0, -R, 1.0};
//+
Circle(1) = {2, 1, 4};
//+
Circle(2) = {4, 1, 3};
//+
Circle(3) = {3, 1, 5};
//+
Circle(4) = {5, 1, 2};
//+
Circle(5) = {6, 1, 8};
//+
Circle(6) = {8, 1, 7};
//+
Circle(7) = {7, 1, 9};
//+
Circle(8) = {9, 1, 6};
//+
Curve Loop(1) = {8, 5, 6, 7};
//+
Curve Loop(2) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1, 2};
//+
Curve Loop(3) = {2, 3, 4, 1};
//+
Extrude {L, 0, 0} {
  Surface{1}; Layers {5};
}
//+
Recursive Delete {
  Point{1}; 
}
//+
Physical Surface("I", 25) = {4, 5, 2, 3};
//+
Physical Surface("U", 26) = {10, 9, 6, 7, 8, 1};
