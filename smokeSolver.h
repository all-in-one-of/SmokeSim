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
    /* Laplacian Matrix */
    Eigen::SparseMatrix<fReal> Laplacian;
    /* grid attribute name : pointer to grid attribute*/
    std::map<std::string, SmokeQuantity*> centerAttr;
    std::map<std::string, SmokeQuantity*> faceAttr;
    std::map<std::string, SmokeQuantity*> attributeTable;
    std::map<std::string, GridStash*> gridStashTable;

    /* returns gridType at specified grid cell */
    gridType getGridTypeAt(size_t x, size_t y, size_t z);
    /* array of gridtypes at cell locations */
    gridType* gridTypes;

    /* advection */
    void advection(fReal dt)
    /* apply forces */
    void force(fReal dt);
    /* pressure projection */
    void precomputeLaplacian();
    void projection(fReal dt);

    /* attribute storage */
    void addCenterAttr(std::string name, fReal xOffset = 0.5, fReal yOffset = 0.5, fReal zOffset = 0.5);
    void addFaceAttr(std::string name, fReal xOffset, fReal yOffset, fReal zOffset);
    void addGridStash(std::string name);
    
    SmokeQuantity* getAttributeNamed(std::string name);
    SmokeQuantity* operator[](std::string name);
    void swapBuffers();
    void swapFaceBuffers();
    void swapCenterBuffers();

public:
    SmokeSolver(size_t Nx, size_t Ny, size_t Nz, fReal h);
    ~SmokeSolver();

    /* advance solver */
    void step(fReal dt);

    /* output bgeo files */
    void write_data_bgeo(const std::string& s, const int frame);
};