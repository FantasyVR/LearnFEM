//
// Created by yupeng on 2021/2/1.
//
#include "Rigidbody.h"
#include "World.h"
void Rigidbody::update()
{
    //TODO: Collision detection

    for(int step = 0; step< numSubSteps; step++)
    {
        //TODO: Semi-Euler update

        //TODO: solve constraints
        for(int ite=0; ite < numIteration; ite++)
        {

        }
        //TODO: update velocities

        //TODO: Post processing of velocities

    }
}
