# include "gridStash.h"

GridStash::GridStash(std::string attrName, size_t Nx, size_t Ny, size_t Nz) :
                                                attrName(attrName), Nx(Nx), Ny(Ny), Nz(Nz)
{
    this->thisStep = new fReal[Nx * Ny * Nz];
}
~GridStash(){
    delete[] thisStep;
}

fReal GridStash::getValueAt(size_t x, size_t y, size_t z)
{
    return this->accessValueAt(x, y, z);
}

fReal& GridStash::accessValueAt(size_t x, size_t y, size_t z)
{
    return this->thisStep[getIndex(x, y, z)];
}

void GridStash::writeValueTo(size_t x, size_t y, size_t z, fReal val)
{
    this->thisStep[getIndex(x, y, z)] = val;
}

size_t GridStash::getIndex(size_t x, size_t y, size_t z)
{
# ifdef DEBUGBUILD
    if (x >= this->Nx || y >= this->Ny || z >= this->Nz)
    {
        std::cerr << "Index out of bound at x: " << x << " y: " << y << " z: " << z << std::endl;
    }
# endif
    return z * (Nx * Ny) + y * (Nx) + x;
}