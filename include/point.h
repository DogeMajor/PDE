#ifndef POINT_H
#define POINT_H
#include <vector>
#include <iostream>
using namespace std;

template <int Dim,typename T> class Point{
public:
    Point();
    Point(vector <T> &val);
    Point(const Point & p);
	const int size() const { return value.size(); }
    vector <T> get_value() const;
    bool operator==(const Point &p) const;
    bool operator!=(const Point &p) const;
    Point <Dim,T>& operator=(const Point &p);
    const T& operator[](int i) const;
	friend const Point<Dim,T> operator+(const Point<Dim,T> &p, const Point<Dim,T> &q) {
		vector<T> loc(p.value.size());
		for (int i = 0; i < p.value.size(); i++) {loc[i] = p.value[i] + q.value[i];}
		return Point(loc);
	}
	friend const Point<Dim,T> operator-(const Point<Dim,T> &p, const Point<Dim,T> &q) {
		vector<T> loc(p.value.size());
		for (int i = 0; i < p.value.size(); i++) {loc[i] = p.value[i] - q.value[i];}
		return Point(loc);
	}
	friend const Point<Dim,T> operator*(const T coeff, const Point<Dim,T> &p) {
		vector<T> loc(p.value.size());
		for (int i = 0; i < p.value.size(); i++) { loc[i] = coeff * p.value[i]; }
		return Point(loc);
	}
    friend const Point <Dim,T> operator*(const Point<Dim,T> &p, const T coeff) {
		return coeff*p;
	}
    void show();
private:
    vector <T> value;

};


template <int Dim, typename T>
Point<Dim,T>::Point() : value(Dim, 0) {}

template <int Dim, typename T>
Point<Dim,T>::Point(vector <T> &val){
    value = val;
}

template <int Dim, typename T>
Point<Dim,T>::Point(const Point &p){
    value = p.value;
}

template <int Dim, typename T>
vector <T> Point<Dim,T>::get_value() const{
    return value;
}

template <int Dim, typename T>
bool Point<Dim,T>::operator==(const Point &p) const{
	return (size() == p.size()) && (value == p.value);
}

template <int Dim, typename T>
bool Point<Dim,T>::operator!=(const Point &p) const{
    return !(*this == p);
}

template <int Dim, typename T>
Point <Dim,T>& Point<Dim,T>::operator=(const Point &p){
    if(p!=*this){value = p.value;}
    return *this;
}

template <int Dim, typename T>
const T& Point<Dim,T>::operator[](int i) const{
    return value[i];
}

template <int Dim, typename T>
void Point<Dim,T>::show(){
    cout << endl;
    for(int i=0; i<value.size(); i++){
        cout << value[i] << " ";
    }
    cout << endl;
}

#endif
