#ifndef WORLD
#define  WORLD
#include "Softbody.h"
class World
{
private:
    Softbody *softbody;
    double dt;
public:
    World() = default;
    World(Softbody* sb):softbody(sb){}
    ~World(){delete softbody;}
    void addSoftbody(Softbody *sb){softbody = sb;}
    void step();
    void setDt(double timestep);
    double getDt();
};

#endif