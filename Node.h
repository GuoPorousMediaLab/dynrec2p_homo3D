#ifndef simplegeometry_node_h
#define simplegeometry_node_h

#include <iostream>

class Layer;

class Node {
public:
    Node(double x, double y, double z1, double z2, int id, Layer* MyLayer) {create_(x, y, z1, z2, id, MyLayer);}  // z1 is below z2
    Node() {create_();}
    double get_x() const {return x_;}
    double get_y() const {return y_;}
    double get_z1() const {return z1_;}
    double get_z2() const {return z2_;}
    int get_ID() const {return ID_;}
    
private:
    double x_;
    double y_;
    double z1_;
    double z2_;
    int ID_;
    Layer* myLayer_;
    void create_(double, double, double, double, int, Layer*);
    void create_();
};
#endif