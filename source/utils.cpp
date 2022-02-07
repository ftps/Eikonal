#include "utils.hpp"

/*
Function that splits the string s into several sub-strings at a character cc
Similar to python's split method
*/
std::vector<std::string> spltString(const std::string& s, const char& cc)
{
    std::vector<std::string> res;
    std::string aux = "";

    for(char c : s){
        if(c == cc){
            if(aux == "") continue;
            res.push_back(aux);
            aux = "";
        }
        else{
            aux += c;
        }
    }
    if(aux != "") res.push_back(aux);

    return res;
}

/*
Function that will read a line of a file fp
If a number of lines n is given, the function
will ignore the first n lines before reading
*/
std::string readLine(std::fstream& fp, const int& n)
{
    std::string res;

    for(int i = 0; i < n; ++i){
        std::getline(fp, res);
    }
    std::getline(fp, res);
    
    return res;
}

/*
Function that performs a trapezoid integraion of a set
of pairXY (std::vector<std::pair<double, double>>) data
*/
double trapInt(const pairXY& data)
{
    double res = 0;
    for(int i = 1; i < (int)data.size(); ++i){
        res += 0.5*(data[i].second + data[i-1].second)*(data[i].first - data[i-1].first);
    }

    return res;
}

/*
Function that returns the closest approach of a line to a point p.
The line is defined by a starting point r and a velocity or
direction v.
*/
Vector3D closeApproach(const Vector3D& p, const Vector3D& r, const Vector3D& v)
{
    return r + v*((v*(p - r))/(v*v));
}

/*
Function that converts a char* to a std::string
*/
std::string convertCharP(const char* s)
{
    uint i = 0;
    while(s[i] != '\0'){
        ++i;
    }

    return std::string(s, i);
}

/*
Function that removes the extension of a string
*/
std::string remExt(const std::string& s)
{
    uint loc = s.find_last_of(".");
    return s.substr(0, loc);
}