# include "smokeQuantity.h"

SmokeQuantity::SmokeQuantity(std::string attrName, size_t Nx, size_t Ny, size_t Nz, fReal h,
                        fReal xOffset, fReal yOffset, fReal zOffset) :
                        attrName(attrName), Nx(Nx), Ny(Ny), Nz(Nz), h(h), invH(1.0 / h),
                        xOffset(xOffset), yOffset(yOffset), zOffset(zOffset)
{
    this->thisStep = new fReal[Nx * Ny * Nz];
    this->nextStep = new fReal[Nx * Ny * Nz];
}

SmokeQuantity::~SmokeQuantity()
{
    delete[] thisStep;
    delete[] nextStep;
}

size_t SmokeQuantity::getNx()
{
    return this->Nx;
}

size_t SmokeQuantity::getNy()
{
    return this->Ny;
}

size_t SmokeQuantity::getNz()
{
    return this->Nz;
}

void SmokeQuantity::swapBuffer()
{
    fReal* tempPtr = this->thisStep;
    this->thisStep = this->nextStep;
    this->nextStep = tempPtr;
}

fReal SmokeQuantity::getValueAt(size_t x, size_t y, size_t z)
{
    return this->accessValueAt(x, y, z);
}

void SmokeQuantity::setValueAt(size_t x, size_t y, size_t z, fReal val)
{
    this->accessValueAt(x, y, z) = val;
}

fReal& SmokeQuantity::accessValueAt(size_t x, size_t y, size_t z)
{
    return this->thisStep[getIndex(x, y, z)];
}

void SmokeQuantity::writeValueTo(size_t x, size_t y, size_t z, fReal val)
{
    this->nextStep[getIndex(x, y, z)] = val;
}

size_t KaminoQuantity::getIndex(size_t x, size_t y, size_t z)
{
# ifdef DEBUGBUILD
    if (x >= this->Nx || y >= this->Ny || z >= this->Nz)
    {
        std::cerr << "Index out of bound at x: " << x << " y: " << y << " z: " << z << std::endl;
    }
# endif
    return z * (Nx * Ny) + y * (Nx) + x;
}

fReal SmokeQuantity::sampleAt(fReal x, fReal y, fReal z)
{
    int xIndex = std::floor(x * invH - this->xOffset);
    int yIndex = std::floor(y * invH - this->yOffset);
    int zIndex = std::floor(z * invH - this->zOffset);

    size_t lowerX = phiIndex < 0 ? this->nPhi - 1 : phiIndex % nPhi;
    size_t upperX = lowerX + 1;
    upperX = upperX >= nPhi ? 0 : upperX;

    size_t lowerY = thetaIndex < 0 ? 0 : thetaIndex;
    size_t upperY = lowerY + 1;

    fReal lowerLeft = getValueAt(lowerX, lowerY);
    fReal upperLeft = getValueAt(lowerX, upperY);
    fReal lowerRight = getValueAt(upperX, lowerY);
    fReal upperRight = getValueAt(upperX, upperY);

    fReal alphaX = x - static_cast<fReal>(std::floor(x));
    fReal alphaY = y - static_cast<fReal>(std::floor(y));

    fReal lerpedLower = KaminoLerp<fReal>(lowerLeft, lowerRight, alphaX);
    fReal lerpedUpper = KaminoLerp<fReal>(upperLeft, upperRight, alphaX);
    fReal lerped = KaminoLerp<fReal>(lerpedLower, lerpedUpper, alphaY);

    return lerped;
}