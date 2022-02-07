#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>
#include <array>
#include <limits>
#include <algorithm>

#define LOG {std::cout << "IN LINE " << __LINE__ << " OF FILE " << __FILE__ << std::endl; fflush(stdout);}

template<typename T, std::size_t N>
using Vector = std::array<T,N>;
using Vector3D = Vector<double, 3>;
using pairXY = std::vector<std::pair<double, double>>;
using uint = unsigned int;


std::vector<std::string> spltString(const std::string& s, const char& cc = ' ');
std::string readLine(std::fstream& fp, const int& = 0);
double trapInt(const pairXY& data);
Vector3D closeApproach(const Vector3D& p, const Vector3D& r, const Vector3D& v);
std::string convertCharP(const char* s);
std::string remExt(const std::string& s);

/*
Auxiliary functions for the algorithms
*/

// print operator for std::pair
template<typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::pair<T,U>& p)
{
    os << "{" << p.first << ", " << p.second << "}";
    return os;
}

// print operator for std::vector
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    os << "[ ";
    for(T e : v){
        os << e << " ";
    }
    os << "]";

    return os;
}

// Function to square a number
template<typename T>
inline T sqr(const T& x)
{
    return x*x;
}

// Function that returns the sign of a value (1 if positive, -1 if negative and 0 if null)
template<typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}


/*
Function to find an element id on a C style array or general container
Returns the index of the element in the vector and if it's
not found, returns -1
*/

template<typename T>
int findIn(const T* v, const T& id, const int& size)
{
    for(int i = 0; i < size; ++i){
        if(v[i] == id) return i;
    }
    return -1;
}

template<typename T, typename U>
int findIn(const T& v, const U& id)
{
    typename T::const_iterator it = std::find(v.begin(), v.end(), id);
    if(it == v.end()){
        return -1;
    }
    return std::distance(v.begin(), it);
}

// functions auxiliary for min and max comparisons

template<typename T>
bool cMin(const T& a, const T& b){ return a > b; }

template<typename T>
bool cMax(const T& a, const T& b){ return a < b; }

/*
Functions for obtaining the miniming and maximum value in a vector.
Returns a pair with a bool stating if the vector has any size, and a
type T with the minimum value in the vector (if no size, set to 0).
*/

template<typename T>
std::pair<bool,T> min(const std::vector<T>& v)
{
    if(v.begin() == v.end()) return {false, 0};

    return {true, *max_element(v.begin(), v.end(), cMin<T>)};
}

template<typename T>
std::pair<bool,T> max(const std::vector<T>& v)
{
    if(v.begin() == v.end()) return {false, 0};

    return {true, *max_element(v.begin(), v.end(), cMax<T>)};
}

/*
Function to remove an element in an unordered vector
*/
template<typename T>
void removeIndex(std::vector<T>& v, const uint i)
{
    if(i >= v.size()) return;
    v[i] = v.back();
    v.pop_back();
}

template<typename T>
void removeItem(const T& e, std::vector<T>& v)
{
    int i = findIn(v, e);
    if(i == -1) return;
    v[i] = v.back();
    v.pop_back();
}



/*
Matrix NxN class
*/
template<typename T, std::size_t N>
class Matrix {
    private:
        // container for the matrix - an array of size N with N arrays of size N
        Vector<Vector<T, N>, N> v;
    public:
        Matrix() {};
        Matrix(const Vector<Vector<T,N>, N>& v) : v(v) {};

        /*
        Matrix operators
        */
        
        Vector<T,N>& operator[](const int& n)
        {
            return v[n];
        }

        Vector<T,N> operator[](const int& n) const
        {
            return v[n];
        }

        Vector<T,N> operator*(const Vector<T,N>& r) const
        {
            Vector<T,N> res = {0};

            for(std::size_t i = 0; i < N; ++i){
                for(std::size_t j = 0; j < N; ++j){
                    res[i] += v[i][j]*r[j];
                }
            }

            return res;
        }

        Matrix<T,N> operator*(const Matrix<T,N>& r) const
        {
            Matrix<T,N> res({0});

            for(std::size_t i = 0; i < N; ++i){
                for(std::size_t j = 0; j < N; ++j){
                    for(std::size_t k = 0; k < N; ++k){
                        res[i][j] += v[i][k]*r[k][j];
                    }
                }
            }

            return res;
        }
};

/*
Remaining matrix operators
*/

template<typename T, std::size_t N>
Matrix<T,N> operator^(Vector<T,N> l, Vector<T,N> r)
{
    Matrix<T,N> m;
    for(std::size_t i = 0; i < N; ++i){
        for(std::size_t j = 0; j < N; ++j){
            m[i][j] = l[i]*r[j];
        }
    }

    return m;
}

template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const Matrix<T,N>& r)
{
    os << "[";
    for(std::size_t i = 0; i < N; ++i){
        os << "[ ";
        for(std::size_t j = 0; j < N; ++j){
            os << r[i][j] << " ";
        }
        os << ((i == (N-1)) ? "]" : "]\n");
    }
    os << "]";

    return os;
}

template<typename T>
Matrix<double,2> inverse(const Matrix<T,2>& M)
{
    Matrix<double, 2> res;
    double det = M[0][0]*M[1][1] - M[0][1]*M[1][0];

    res[0][0] = M[1][1]/det;
    res[0][1] = -M[0][1]/det;
    res[1][0] = -M[1][0]/det;
    res[1][1] = M[0][0]/det;

    return res;
}

template<typename T>
Matrix<double,3> inverse(const Matrix<T,3>& M)
{
    Matrix<double, 3> res;
    double idet = 1/(M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
		- M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
		+ M[0][2]*(M[1][0]*M[2][1] - M[2][0]*M[1][1]));

    res[0][0] = (M[1][1] * M[2][2] - M[2][1] * M[1][2]) * idet;
	res[0][1] = (M[0][2] * M[2][1] - M[0][1] * M[2][2]) * idet;
	res[0][2] = (M[0][1] * M[1][2] - M[0][2] * M[1][1]) * idet;
	res[1][0] = (M[1][2] * M[2][0] - M[1][0] * M[2][2]) * idet;
	res[1][1] = (M[0][0] * M[2][2] - M[0][2] * M[2][0]) * idet;
	res[1][2] = (M[1][0] * M[0][2] - M[0][0] * M[1][2]) * idet;
	res[2][0] = (M[1][0] * M[2][1] - M[2][0] * M[1][1]) * idet;
	res[2][1] = (M[2][0] * M[0][1] - M[0][0] * M[2][1]) * idet;
	res[2][2] = (M[0][0] * M[1][1] - M[1][0] * M[0][1]) * idet;

    return res;
}



/*
Vector operators
*/

template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const Vector<T,N>& r)
{
    os << "[ ";
    for(std::size_t i = 0; i < N; ++i){
        os << r[i] << " ";
    }
    os << "]";

    return os;
}

template<typename T, std::size_t N>
Vector<T,N> operator+(const Vector<T,N>& l, const Vector<T,N>& r)
{
    Vector<T,N> res;
    for(std::size_t i = 0; i < N; ++i){
        res[i] = l[i] + r[i];
    }

    return res;
}

template<typename T, std::size_t N>
Vector<T,N> operator+=(Vector<T,N>& l, const Vector<T,N>& r)
{
    for(std::size_t i = 0; i < N; ++i){
        l[i] += r[i];
    }

    return l;
}

template<typename T, std::size_t N>
Vector<T,N> operator-(const Vector<T,N>& l, const Vector<T,N>& r)
{
    Vector<T,N> res;
    for(std::size_t i = 0; i < N; ++i){
        res[i] = l[i] - r[i];
    }

    return res;
}

template<typename T, std::size_t N>
Vector<T,N> operator-(const Vector<T,N>& r)
{
    Vector<T,N> res;
    for(std::size_t i = 0; i < N; ++i){
        res[i] = -r[i];
    }

    return res;
}

template<typename T, std::size_t N>
Vector<T,N> operator-=(Vector<T,N>& l, const Vector<T,N>& r)
{
    for(std::size_t i = 0; i < N; ++i){
        l[i] -= r[i];
    }

    return l;
}

template<typename T, std::size_t N>
T operator*(const Vector<T,N>& l, const Vector<T,N>& r)
{
    T res = 0;
    for(std::size_t i = 0; i < N; ++i){
        res += l[i]*r[i];
    }

    return res;
}

template<typename U, typename T, std::size_t N>
Vector<T,N> operator*(const Vector<T,N>& l, const U& r){
    Vector<T,N> res;
    for(std::size_t i = 0; i < N; ++i){
        res[i] = r*l[i];
    }
    return res;
}

template<typename U, typename T, std::size_t N>
Vector<T,N> operator*(const U& r, const Vector<T,N>& l){
    return l*r;
}

template<typename U, typename T, std::size_t N>
Vector<T,N> operator/(const Vector<T,N>& l, const U& r){
    Vector<T,N> res;
    for(std::size_t i = 0; i < N; ++i){
        res[i] = l[i]/r;
    }
    return res;
}

/*
2-norm of a vector
*/
template<typename T, std::size_t N>
T modV(const Vector<T,N>& l)
{
    return sqrt(l*l);
}

/*
Cross products between two 3D vectors
*/
template<typename T>
Vector<T,3> cross(const Vector<T,3>& l, const Vector<T,3>& r)
{
    Vector<T,3> res;
    for(int i = 0; i < 3; ++i){
        res[i] = l[(i+1)%3]*r[(i+2)%3] - l[(i+2)%3]*r[(i+1)%3];
    }

    return res;
}

#endif