# pragma once

# include "globalIncludes.h"
# include "smokeQuantity.h"
# include "gridStash.h"

class SmokeSolver
{
private:
    /* grid dimensions */
    size_t Nx, Ny, Nz;
    /* grid spacing */
    fReal h;
    fReal invH;
    /* Laplacian Matrix */
    Eigen::SparseMatrix<fReal> Laplacian;
    /* preconditioner for incomplete Cholesky */
    Eigen::VectorXd precon;
    /* force parameters */
    fReal alpha;
    fReal beta;
    fReal epsilon;
    /* velocity source initialization speed */
    fReal sourceSpeed;
    
    /* grid attribute name : pointer to grid attribute*/
    std::map<std::string, SmokeQuantity*> centerAttr;
    std::map<std::string, SmokeQuantity*> faceAttr;
    std::map<std::string, SmokeQuantity*> attributeTable;
    std::map<std::string, GridStash*> gridStashTable;

    /* returns gridType at specified grid cell */
    gridType getGridTypeAt(size_t x, size_t y, size_t z);
    gridType getGridTypeAt(fReal x, fReal y, fReal z);
    /* array of gridtypes at cell locations */
    gridType* gridTypes;
    /* get index */
    inline size_t getIndex(size_t x, size_t y, size_t z);

    /* advection */
    void advection(fReal dt);

    /* apply forces */
    void force(fReal dt);
    void buoyancyForce();
    void vorticity();
    void forceVelocityUpdate(fReal dt);

    /* pressure projection */
    void precomputeLaplacian();
    void testLaplacian();
    void projection(fReal dt);
    void testBVector(Eigen::VectorXd& b);
    void fillDivergence(Eigen::VectorXd& b, fReal uSolid, fReal vSolid, fReal wSolid);
    void testDivergenceFree(Eigen::VectorXd& b);
    void updateVelWithPressure(SmokeQuantity* speed, GridStash* p, fReal scaleP);
    void PCG(Eigen::VectorXd& p, Eigen::VectorXd& b);
    void pressureSolve(Eigen::VectorXd& p, Eigen::VectorXd& b);
    void computePreconditioner();
    void applyPreconditioner(Eigen::VectorXd& z, Eigen::VectorXd& r);

    /* attribute storage */
    void addCenterAttr(std::string name, fReal xOffset, fReal yOffset, fReal zOffset);
    void addFaceAttr(std::string name, fReal xOffset, fReal yOffset, fReal zOffset);
    void addCenterGridStash(std::string name);
    void addFaceGridStash(std::string name, size_t xO, size_t yO, size_t zO);
    
    /* convenient access to maps stored in solver */
    SmokeQuantity* operator[](std::string name);

    /* for buffer swapping */
    void swapBuffers();
    void swapFaceBuffers();
    void swapCenterBuffers();

    /* for averaging quantities across faces / centers */
    void faceToCenterVelocity(GridStash* centerStash, SmokeQuantity* vel, size_t xO, size_t yO, size_t zO);
    void centerToFaceStash(GridStash* faceStash, GridStash* centerStash, size_t xO, size_t yO, size_t zO);

    /* boundary initialization */
    void initializeBoundary();

    /* set all sources */
    void setSources();
    /* reset source velocity */
    void setSourceVelocity(SmokeQuantity* speed, size_t x, size_t y, size_t z, fReal val);
    /* reset source temperature */
    void setSourceTemperature(size_t x, size_t y, size_t z, fReal val);
    /* reset source density */
    void setSourceDensity(size_t x, size_t y, size_t z, fReal val);

public:
    SmokeSolver(size_t Nx, size_t Ny, size_t Nz, fReal h, fReal alpha, fReal beta, fReal epsilon, fReal sourceSpeed);
    ~SmokeSolver();

    /* advance solver */
    void step(fReal dt);
    /* get pointer to MAC grid arrays */
    SmokeQuantity* getAttributeNamed(std::string name);
    /* get pointer to stashed arrays */
    GridStash* getStashNamed(std::string name);

    /* output bgeo files */
    void write_data_bgeo(const std::string& s, const int frame);
};