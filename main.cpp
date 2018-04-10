# include "globalIncludes.h"
# include "smokeSolver.h"
# include "smokeQuantity.h"

/* grid specs */
const size_t Nx = 50;
const size_t Ny = 50;
const size_t Nz = 50;
fReal h = 0.1;

/* time step and frame rate */
const float dt = 0.005;
const float DT = 1.0 / 24.0;
const size_t frames = 50;

/* file output destination */
const std::string filepath = "output/frame";

int main(int argc, char** argv)
{
    SmokeSolver solver();
    float T = 0.0;              // simulation time
    for(int i = 1; i <= frames; i++){
        while(T < i*DT){
            //solver.stepForward(dt);
            T += dt;
        }
        //solver.stepForward(dt + i*DT - T);
        T = i*DT;
        //solver.write_data_bgeo(filepath, i);
    }
    return 0;
}
