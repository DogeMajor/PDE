#ifndef POINT_H
#define POINT_H
#include <vector>
#include <iostream>
using namespace std;

template <typename T> class Point{
public:
    Point();
    Point(vector <T> val);
    Point(Point<T> & p);
    int get_dimension() const;
    vector <T> get_value() const;
    bool operator==(const Point<T> &p) const;
    bool operator!=(const Point<T> &p) const;
    Point <T> operator=(const Point<T> &p);
    T operator[](int i);
    void show();
private:
    vector <T> value;

};


template <typename T>
Point<T>::Point(){
    //value = {0};
}

template <typename T>
Point<T>::Point(vector <T> val){
    value = val;
}

template <typename T>
Point<T>::Point(Point <T> &p){
    value = p.value;
}

template <typename T>
int Point<T>::get_dimension() const{
    return value.size();
}

template <typename T>
vector <T> Point<T>::get_value() const{
    return value;
}

template <typename T>
bool Point<T>::operator==(const Point<T> &p) const{
    return (value == p.value);
}

template <typename T>
bool Point<T>::operator!=(const Point<T> &p) const{
    return !(*this == p);
}

template <typename T>
Point <T> Point<T>::operator=(const Point<T> &p){
    if(p!=*this) value == p.value;
    return *this;
}

template <typename T>
T Point<T>::operator[](int i){
    return value[i];
}

template <typename T>
void Point<T>::show(){
    cout << endl;
    for(int i=0; i<value.size(); i++){
        cout << value[i] << " ";
    }
    cout << endl;
}

#endif
