# include "globalIncludes.h"
# include "smokeQuantity.h"
# include "smokeParticles.h"
# include "smokeSolver.h"

/* grid specs */
const size_t Nx = 20;
const size_t Ny = 20;
const size_t Nz = 20;
/* for a 1x1x1 cube */
fReal h = 1.0 / Nx;

/* time step and frame rate */
/* CFL condition: dt <= 5*dx / u_max */
/* dt < DT */
const float dt = 0.005;
const float DT = 1.0 / 24.0;
const size_t frames = 250;

/* particle birth rate */
const size_t particleRate = 10;

/* source location */
const size_t Sx = Nx / 2;
const size_t Sy = Ny / 5;
const size_t Sz = Nz / 2;

/* file output destination */
const std::string filepath = "output/frame";

int main(int argc, char** argv)
{
    SmokeSolver solver(Nx, Ny, Nz, h);
    SmokeParticles particles(h);
    particles.addParticles(Sx, Sy, Sz, particleRate);
    particles.write_data_bgeo(filepath, 1);

    SmokeQuantity* u = solver.getAttributeNamed("u");
    SmokeQuantity* v = solver.getAttributeNamed("v");
    SmokeQuantity* w = solver.getAttributeNamed("w");

    float T = 0.0;              // simulation time
    
    std::cout << "writing frame 0" << std::endl;

    for(int i = 1; i <= frames; i++){
        while(T < i*DT){
            solver.setSourceVelocity(Sx, Sy, Sz, 2.0);
            solver.setSourceTemperature(Sx, Sy, Sz, 1.0);
            solver.setSourceDensity(Sx, Sy, Sz, 1.0);
            solver.step(dt);
            particles.updatePositions(u, v, w, dt);
            particles.addParticles(Sx, Sy, Sz, particleRate);
            T += dt;
        }
        solver.setSourceVelocity(Sx, Sy, Sz, 2.0);
        solver.setSourceTemperature(Sx, Sy, Sz, 1.0);
        solver.setSourceDensity(Sx, Sy, Sz, 1.0);
        solver.step(dt + i*DT - T);
        particles.updatePositions(u, v, w, dt);
        particles.addParticles(Sx, Sy, Sz, particleRate);
        
        T = i*DT;
        
        std::cout << "writing frame " << i << std::endl;
        particles.write_data_bgeo(filepath, i);
    }
    return 0;
}
