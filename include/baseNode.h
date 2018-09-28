#ifndef BASENODE_H
#define BASENODE_H
#include <iostream>


class BaseNode{

public:
    BaseNode(){}
    //BaseNode(const BaseNode &a);//copy constructor
    virtual ~BaseNode(){}
    //virtual void set_index(int ind) = 0;
    //virtual void set_shared_elements(int shared_els) = 0;
    //virtual int get_index() const = 0;
    //virtual int get_shared_elements() const = 0;
    //virtual int get_BaseNode_amount() const = 0;
    //BaseNode& operator=(const BaseNode &a);
    //bool operator== (const BaseNode &a) const;
    //bool operator!=(const BaseNode &a) const;
    virtual void show() const = 0;

};

#endif

