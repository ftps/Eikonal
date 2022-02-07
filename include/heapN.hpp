#ifndef HEAP_HPP
#define HEAP_HPP

#include <functional>
#include "utils.hpp"

/*
Function meant to represent the weight of a heap element
with the ID as input parameter
*/
using functionW = std::function<double(const int&)>;

/*
Data structure that implements a priority queue using a binary tree
Implementation in standard C++17, faster performance than C-style
*/
class HeapN {
    public:
        /*
        Heap parameters:
        - w: function takes in an int (heap id) and return the weight of that element
        - v: vector with heap ordered elements
        */
        const functionW w;
        std::vector<int> v;

    public:
        HeapN(const functionW& w);
        void reserve(const size_t& s);
        bool isEmpty();
        uint size();
        void addNew(const int& n);
        int extractRoot();
        void siftUp(const int& i_f = -1);
        void siftDown(const int& i_f = -1);
        void printHeap();
};

#endif