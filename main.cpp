# include "globalIncludes.h"
# include "smokeQuantity.h"
# include "smokeParticles.h"
# include "smokeSolver.h"

/* grid specs */
// recommend odd numbers for centrally placed sources
const size_t Nx = 25;
const size_t Ny = 25;
const size_t Nz = 25;
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
const fReal epsilon = 0.02;

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

    size_t Sx1 = Nx / 2;
    size_t Sy1 = 2;
    size_t Sz1 = Nz / 2;

    size_t Sx2 = 3 * Nx / 5;
    size_t Sy2 = 2;
    size_t Sz2 = 3 * Nz / 5;

    particles.addParticles(Sx1, Sy1, Sz1, particleRate);
    // particles2.addParticles(Sx2, Sy2, Sz2, particleRate);

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
            particles.addParticles(Sx1, Sy1, Sz1, particleRate);
            // particles2.updateParticles(u, v, w, dt);
            // particles2.addParticles(Sx2, Sy2, Sz2, particleRate);

            T += dt;
        }
        solver.step(dt + i*DT - T);
        particles.updateParticles(u, v, w, dt);
        particles.addParticles(Sx1, Sy1, Sz1, particleRate);
        // particles2.updateParticles(u, v, w, dt);
        // particles2.addParticles(Sx2, Sy2, Sz2, particleRate);

        T = i*DT;
        
        std::cout << "writing frame " << i << std::endl;
# ifdef PARTIO

        particles.write_data_bgeo(filepathP, i);
        // particles2.write_data_bgeo(filepathP2, i);

        // solver.write_data_bgeo(filepathG, i);
# endif
    }
    return 0;
}