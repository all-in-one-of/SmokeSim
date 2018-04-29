# include "globalIncludes.h"
# include "smokeQuantity.h"
# include "smokeParticles.h"
# include "smokeSolver.h"

/* grid specs */
const size_t Nx = 25;
const size_t Ny = 25;
const size_t Nz = 25;
/* for a 1x1x1 cube */
fReal h = 1.0 / Nx;

/* time step and frame rate */
/* CFL condition: dt <= 5*dx / u_max */
/* dt < DT */
const float dt = 0.005;
const float DT = 1.0 / 24.0;
const size_t frames = 100;

/* particle birth rate */
const size_t particleRate = 15;

// NOTE: Sources defined in SmokeSolver > setSources()

/* file output destination */
const std::string filepath = "output/frame";

int main(int argc, char** argv)
{
    size_t Sx = Nx / 2;
    size_t Sy = 1;
    size_t Sz = Nz / 2;

    SmokeSolver solver(Nx, Ny, Nz, h);
    SmokeParticles particles(h);
    particles.addParticles(Sx, Sy, Sz, particleRate);
    particles.write_data_bgeo(filepath, 1);

    SmokeQuantity* u = solver.getAttributeNamed("u");
    SmokeQuantity* v = solver.getAttributeNamed("v");
    SmokeQuantity* w = solver.getAttributeNamed("w");

    float T = 0.0;              // simulation time
    
    std::cout << "writing frame 1" << std::endl;

    for(int i = 2; i <= frames; i++){
        while(T < i*DT){
            solver.step(dt);
            particles.updateParticles(u, v, w, dt);
            particles.addParticles(Sx, Sy, Sz, particleRate);
            T += dt;
        }
        solver.step(dt + i*DT - T);
        particles.updateParticles(u, v, w, dt);
        particles.addParticles(Sx, Sy, Sz, particleRate);
        
        T = i*DT;
        
        std::cout << "writing frame " << i << std::endl;
        particles.write_data_bgeo(filepath, i);
    }
    return 0;
}
