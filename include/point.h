#ifndef POINT_H
#define POINT_H
#include <vector>
#include <iostream>
using namespace std;

template <class T> class Point{
public:
    Point(vector <T> val);
    void set_index(int ind);
    int get_index() const;
    vector <T> get_value() const;
private:
    vector <T> value;
    int index;

};


template <class T>
Point<T>::Point(vector <T> val){
    value = val;
    index = -1;
}

template <class T>
void Point<T>::set_index(int ind){
    index = ind;
}

template <class T>
int Point<T>::get_index() const{
    return index;
}

template <class T>
vector <T> Point<T>::get_value() const{
    return value;
}


#endif
