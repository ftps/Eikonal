#include "heapN.hpp"

/*
Heap constructor, takes in the max size and the weight function
Allocates memory for the heap vector
*/
HeapN::HeapN(const functionW& w) : w(w), v(0) { }

/*
Reserves memory in the heap for s elements
Useful when a theoretical maximum is known
to avoid unecessary allocation of memory
*/
void HeapN::reserve(const size_t& s)
{
    v.reserve(s);
}

// returns true only if the heap is empty
bool HeapN::isEmpty()
{
    return !v.size();
}

// returns current heap size
uint HeapN::size()
{
    return v.size();
}

/*
Adds a new element of id n to the heap
*/
void HeapN::addNew(const int& n)
{
    v.emplace_back(n);
    siftUp();
}

/*
Extracts the root of the heap, or the element with smallest weight
When this is done, the last element of the heap is dragged to the
root and sifted down towards its correct position
*/
int HeapN::extractRoot()
{
    int root = v[0];

    v[0] = v.back();
    v.pop_back();
    siftDown();

    return root;
}

/*
Heap algorithm to place a newly inserted element in its place,
or one that has had its weight reduced
It compares the weight of the current selected element with the one
before it and proceeds to swtich them until either finding an element
that is indeed smaller or until reaching the root (smallest element)
Takes in the id of the element in the heap, or if none is given it will
use the last element in the heap (used in insertion)
*/
void HeapN::siftUp(const int& i_f)
{
    int i, j, iaux;

    if(i_f == -1) i = v.size()-1;
    else i = findIn(v, i_f);

    while(i > 0){
        j = (i-1)/2;
        if(w(v[i]) > w(v[j])) break;

        iaux = v[i];
        v[i] = v[j];
        v[j] = iaux;
        i = j;
    }
}

/*
Heap algorithm that does the reverse of the siftUp algorithm
Compares the weight of the current selected element with the
weight of the next two elements and performs a switch with the
smallest one if there is one, or until reaching the end of the
heap. Starts at the given id number or, if none is given, at
the root of the heap (used in removal)
*/
void HeapN::siftDown(const int& i_f)
{
    int i, j, j1, j2, iaux;

    if(i_f == -1) i = 0;
    else i = findIn(v, i_f);

    while(i < (int)v.size()){
        j1 = 2*(i+1) - 1;
        if(j1 >= (int)v.size()) break;
        j2 = 2*(i+1);
        if(j2 >= (int)v.size()) j = j1;
        else j = (w(v[j1]) < w(v[j2])) ? j1 : j2;
        
        if(w(v[i]) < w(v[j])) break;

        iaux = v[i];
        v[i] = v[j];
        v[j] = iaux;
        i = j;
    }
}

/*
Prints the heap vector with the ordered id
and weight side by side
*/
void HeapN::printHeap()
{
    std::cout << "Heap: " << std::endl;
    for(const int& i : v){
        std::cout << "(" << i << ", " << w(i) << ")" << ((i == v.back()) ? "\n" : ", ");
    }
}