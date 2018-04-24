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

    /* pressure projection */
    void precomputeLaplacian();
    void testLaplacian();
    void projection(fReal dt);
    void testBVector(Eigen::VectorXd& b);
    void fillDivergence(Eigen::VectorXd& b, fReal uSolid, fReal vSolid, fReal wSolid);
    void testDivergenceFree(Eigen::VectorXd& b);

    /* attribute storage */
    void addCenterAttr(std::string name, fReal xOffset = 0.5, fReal yOffset = 0.5, fReal zOffset = 0.5);
    void addFaceAttr(std::string name, fReal xOffset, fReal yOffset, fReal zOffset);
    void addCenterGridStash(std::string name);
    void addFaceGridStash(std::string name, int xO, int yO, int zO);
    
    /* convenient access to maps stored in solver */
    SmokeQuantity* operator[](std::string name);

    /* for mass buffer swapping */
    void swapBuffers();
    void swapFaceBuffers();
    void swapCenterBuffers();

    /* for averaging quantities across faces / centers */
    void averageVelocity(GridStash* avgStash, SmokeQuantity* attr, int xO, int yO, int zO);
    void averageCenterStash(GridStash* avgStash, GridStash* stash, int xO, int yO, int zO);

    /* boundary initialization */
    void initializeBoundary();

public:
    SmokeSolver(size_t Nx, size_t Ny, size_t Nz, fReal h);
    ~SmokeSolver();

    /* advance solver */
    void step(fReal dt);
    /* reset source velocity */
    void setSourceVelocity(size_t x, size_t y, size_t z, fReal val);
    /* reset source temperature */
    void setSourceTemperature(size_t x, size_t y, size_t z, fReal val);
    /* reset source density */
    void setSourceDensity(size_t x, size_t y, size_t z, fReal val);
    /* get pointer to MAC grid arrays */
    SmokeQuantity* getAttributeNamed(std::string name);
    /* get pointer to stashed arrays */
    GridStash* getStashNamed(std::string name);

    /* output bgeo files */
    void write_data_bgeo(const std::string& s, const int frame);
};