#include "Layer.h"
#include "Node.h"

void Node::create_(double x, double y, double z1, double z2, int id, Layer* myLayer) {
    x_ = x;
    y_ = y;
    z1_ = z1;
    z2_ = z2;
    ID_ = id;
    myLayer_ = myLayer;
    
}
void Node::create_() {
    x_ = 0.0;
    y_ = 0.0;
    z1_ = 0.0;
    z2_ = 0.0;
    ID_ = 0;
}

