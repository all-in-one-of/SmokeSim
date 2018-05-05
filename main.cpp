# include "globalIncludes.h"
# include "smokeQuantity.h"
# include "smokeParticles.h"
# include "smokeSolver.h"

/* grid specs */
// recommend odd numbers for centrally placed sources
const size_t Nx = 21;
const size_t Ny = 21;
const size_t Nz = 21;
/* for a 1x1x1 cube */
const fReal h = 1.0 / static_cast <fReal> (Nx);

/* time step and frame rate */
/* CFL condition: dt <= 5*dx / u_max */
/* dt < DT */
const float dt = 0.005;
const float DT = 1.0 / 24.0;
const size_t frames = 200;

/* particle birth rate */
const size_t particleRate = 200;

/* source speed */
const fReal sourceSpeed = 0.75;

/* gravity */
const fReal alpha = 0.08;
/* buoyancy */
const fReal beta = 0.37;
//const fReal beta = 0.17;
/* vorticity gain */
const fReal epsilon = 0.03;

// NOTE: Sources defined in SmokeSolver > setSources()

/* file output destination */
const std::string filepathP = "particles/frame";
const std::string filepathP2 = "particles2/frame";
const std::string filepathP3 = "particles3/frame";
const std::string filepathG = "grid/frame";

int main(int argc, char** argv)
{
    SmokeSolver solver(Nx, Ny, Nz, h, alpha, beta, epsilon, sourceSpeed);

    SmokeParticles particles(h);
    // SmokeParticles particles2(h);

    particles.addParticles(Nx / 2, 2, Nz / 2, particleRate);
    // particles2.addParticles(3 * Nx / 5, 2, 3 * Nz / 5, particleRate);

# ifdef PARTIO
    particles.write_data_bgeo(filepathP, 1);
    // particles2.write_data_bgeo(filepathP2, 1);

    solver.write_data_bgeo(filepathG, 1);
# endif

    SmokeQuantity* u = solver.getAttributeNamed("u");
    SmokeQuantity* v = solver.getAttributeNamed("v");
    SmokeQuantity* w = solver.getAttributeNamed("w");

    float T = 0.0;              // simulation time
    
    std::cout << "writing frame 1" << std::endl;

    for(int i = 2; i <= frames; i++){
        while(T < i*DT){
            solver.step(dt);
            particles.updateParticles(u, v, w, dt);
            particles.addParticles(Nx / 2, 2, Nz / 2, particleRate);
            // particles2.updateParticles(u, v, w, dt);
            // particles2.addParticles(3 * Nx / 5, 2, 3 * Nz / 5, particleRate);

            T += dt;
        }
        solver.step(dt + i*DT - T);
        particles.updateParticles(u, v, w, dt);
        particles.addParticles(Nx / 2, 2, Nz / 2, particleRate);
        // particles2.updateParticles(u, v, w, dt);
        // particles2.addParticles(3 * Nx / 5, 2, 3 * Nz / 5, particleRate);

        T = i*DT;
        
        std::cout << "writing frame " << i << std::endl;
# ifdef PARTIO

        particles.write_data_bgeo(filepathP, i);
        // particles2.write_data_bgeo(filepathP2, i);

        solver.write_data_bgeo(filepathG, i);
# endif
    }
    return 0;
}