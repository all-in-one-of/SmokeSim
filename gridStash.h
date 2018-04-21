# pragma once

# include "globalIncludes.h"

class GridStash
{
private:
    /* attribute name */
    std::string attrName;
    /* grid dimensions */
    size_t Nx, Ny, Nz;

    /* 3D storage array */
    fReal* thisStep;

    /* get index of buffer array at grid cell location */
    size_t getIndex(size_t x, size_t y, size_t z);

public:
    GridStash(std::string attrName, size_t Nx, size_t Ny, size_t Nz);
    ~GridStash();

    /* get read attribute at grid cell */
    fReal getValueAt(size_t x, size_t y, size_t z);
    /* write attribute value to THIS buffer */
    void writeValueTo(size_t x, size_t y, size_t z, fReal val);
    /* get write attribute at grid cell */
    fReal& accessValueAt(size_t x, size_t y, size_t z);
};