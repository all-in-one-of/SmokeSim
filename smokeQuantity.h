# pragma once

# include "globalIncludes.h"
# define DEBUGBUILD

class SmokeQuantity
{
private:
    /* attribute name */
    std::string attrName;

    /* grid dimensions */
    size_t Nx, Ny, Nz;
    /* grid size */
    fReal h;
    fReal invH;

    /* position offset */
    fReal xOffset, yOffset, zOffset;

    /* get index */
    inline size_t getIndex(size_t x, size_t y, size_t z);

    /* double buffer */
    fReal* thisStep;
    fReal* nextStep;

public:
    SmokeQuantity(std::string attrName, size_t Nx, size_t Ny, size_t Nz, fReal h,
                fReal xOffset, fReal yOffset, fReal zOffset);
    ~SmokeQuantity();

    /* buffer swap */
    void swapBuffer();

    /* get grid size values */
    size_t getNx();
    size_t getNy();
    size_t getNz();

    /* get coordinates of grid locations */
    fReal getXCoordAtIndex(size_t x);
    fReal getYCoordAtIndex(size_t y);
    fReal getZCoordAtIndex(size_t z);

    /* get read attribute at grid cell */
    fReal getValueAt(size_t x, size_t y, size_t z);
    /* set attribute at grid cell for initialization */
    void setValueAt(size_t x, size_t y, size_t z, fReal val);
    /* write attribute value to NEXT buffer */
    void writeValueTo(size_t x, size_t y, size_t z, fReal val);
    /* get write attribute at grid cell */
    fReal& accessValueAt(size_t x, size_t y, size_t z);
    /* interpolation */
    fReal sampleAt(fReal x, fReal y, fReal z);
    /* get name of this quantity */
    std::string getName();
};