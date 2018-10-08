#ifndef POINT_H
#define POINT_H
#include <vector>
#include <iostream>
using namespace std;

template <typename T> class Point{
public:
    Point();
    Point(vector <T> &val);
    Point(const Point & p);
    int get_dimension() const;
    vector <T> get_value() const;
    bool operator==(const Point &p) const;
    bool operator!=(const Point &p) const;
    Point <T>& operator=(const Point &p);
    T operator[](int i);
    //const T& operator*(T coeff);
    //friend Point<T> operator*(T coeff, const Point p){
    //        return p*coeff;//Has to be defined inside the class brackets!!!
    //}
    //vector <T> operator()();replace get_value perhaps
    void show();
private:
    vector <T> value;

};


template <typename T>
Point<T>::Point(){
    //value = {0};
}

template <typename T>
Point<T>::Point(vector <T> &val){
    value = val;
}

template <typename T>
Point<T>::Point(const Point &p){
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
bool Point<T>::operator==(const Point &p) const{
    return (value == p.value);
}

template <typename T>
bool Point<T>::operator!=(const Point &p) const{
    return !(*this == p);
}

template <typename T>
Point <T>& Point<T>::operator=(const Point &p){
    if(p!=*this){
        //value.resize(p.size());
        value = p.value;
    }
    return *this;
}

template <typename T>
T Point<T>::operator[](int i){
    return value[i];
}
/*
template <typename T>
const T& Point<T>::operator*(T coeff){
    for(int i=0; i<value.size(); i++){value[i] = coeff*value[i];}
    return *this;
}*/

template <typename T>
void Point<T>::show(){
    cout << endl;
    for(int i=0; i<value.size(); i++){
        cout << value[i] << " ";
    }
    cout << endl;
}

#endif
