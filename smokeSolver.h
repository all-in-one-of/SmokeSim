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
    void projection(fReal dt);

    /* attribute storage */
    void addCenterAttr(std::string name, fReal xOffset = 0.5, fReal yOffset = 0.5, fReal zOffset = 0.5);
    void addFaceAttr(std::string name, fReal xOffset, fReal yOffset, fReal zOffset);
    void addGridStash(std::string name, int offset);
    
    /* convenient access to maps stored in solver */
    SmokeQuantity* operator[](std::string name);

    /* for mass buffer swapping */
    void swapBuffers();
    void swapFaceBuffers();
    void swapCenterBuffers();

    /* for averaging quantities across faces / centers */
    void averageAttribute(GridStash* avgStash, SmokeQuantity* attr);

    /* boundary initialization */
    void initializeBoundary();

public:
    SmokeSolver(size_t Nx, size_t Ny, size_t Nz, fReal h);
    ~SmokeSolver();

    /* advance solver */
    void step(fReal dt);
    /* reset source velocity */
    void setSourceVelocity(size_t x, size_t y, size_t z, fReal val);
    /* get pointer to MAC grid arrays */
    SmokeQuantity* getAttributeNamed(std::string name);

    /* output bgeo files */
    void write_data_bgeo(const std::string& s, const int frame);
};