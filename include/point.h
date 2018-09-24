#ifndef POINT_H
#define POINT_H
#include <vector>
#include <iostream>
using namespace std;

template <typename T> class Point{
public:
    Point(vector <T> val);
    void set_index(int ind);
    int get_index() const;
    vector <T> get_value() const;
    void show();
private:
    vector <T> value;
    int index;

};


template <typename T>
Point<T>::Point(vector <T> val){
    value = val;
    index = -1;
}

template <typename T>
void Point<T>::set_index(int ind){
    index = ind;
}

template <typename T>
int Point<T>::get_index() const{
    return index;
}

template <typename T>
vector <T> Point<T>::get_value() const{
    return value;
}

template <typename T>
void Point<T>::show(){
    cout << "index:"  << index << endl;
}

#endif
