# include "smokeQuantity.h"

SmokeQuantity::SmokeQuantity(std::string attrName, size_t Nx, size_t Ny, size_t Nz, fReal h,
                        fReal xOffset, fReal yOffset, fReal zOffset) :
                        attrName(attrName), Nx(Nx), Ny(Ny), Nz(Nz), h(h), invH(1.0 / h),
                        xOffset(xOffset), yOffset(yOffset), zOffset(zOffset)
{
    this->thisStep = new fReal[Nx * Ny * Nz]();
    this->nextStep = new fReal[Nx * Ny * Nz]();
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

fReal SmokeQuantity::getXCoordAtIndex(size_t x)
{
    fReal xReal = static_cast<fReal>(x) + xOffset;
    return xReal * this->h;
}

fReal SmokeQuantity::getYCoordAtIndex(size_t y)
{
    fReal yReal = static_cast<fReal>(y) + yOffset;
    return yReal * this->h;
}

fReal SmokeQuantity::getZCoordAtIndex(size_t z)
{
    fReal zReal = static_cast<fReal>(z) + zOffset;
    return zReal * this->h;
}

void SmokeQuantity::swapBuffer()
{
    fReal* tempPtr = this->thisStep;
    this->thisStep = this->nextStep;
    this->nextStep = tempPtr;
}

fReal SmokeQuantity::getValueAt(size_t x, size_t y, size_t z)
{
    if(x >= this->Nx || y >= this->Ny || z >= this->Nz){
        return 0.0;
    }
    else{
        return this->accessValueAt(x, y, z);
    }
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

size_t SmokeQuantity::getIndex(size_t x, size_t y, size_t z)
{
# ifdef DEBUGBUILD
    if (x >= this->Nx || y >= this->Ny || z >= this->Nz)
    {
        std::cerr << "Index out of bound at x: " << x << " y: " << y << " z: " << z << std::endl;
    }
# endif
    return (z * (this->Nx * this->Ny)) + (y * this->Nx) + x;
}

std::string SmokeQuantity::getName()
{
    return this->attrName;
}

// Linear interpolation
// TODO cubic interpolation
fReal Lerp(const fReal &fromEndPoint, const fReal &toEndPoint, fReal ratio)
{
    return (1.0 - ratio) * fromEndPoint + (ratio * toEndPoint);
}

fReal SmokeQuantity::sampleAt(fReal x, fReal y, fReal z)
{
    int xIndex = std::floor(x * this->invH - this->xOffset);
    int yIndex = std::floor(y * this->invH - this->yOffset);
    int zIndex = std::floor(z * this->invH - this->zOffset);

    size_t lowerX = xIndex;// < 0 ? 0 : xIndex;
    size_t upperX = lowerX + 1;
    size_t lowerY = yIndex;// < 0 ? 0 : yIndex;
    size_t upperY = lowerY + 1;
    size_t lowerZ = zIndex;// < 0 ? 0 : zIndex;
    size_t upperZ = lowerZ + 1;

    fReal LHL = getValueAt(lowerX, upperY, lowerZ);
    fReal LLL = getValueAt(lowerX, lowerY, lowerZ);
    fReal LHH = getValueAt(lowerX, upperY, upperZ);
    fReal LLH = getValueAt(lowerX, lowerY, upperZ);
    fReal HHL = getValueAt(upperX, upperY, lowerZ);
    fReal HLL = getValueAt(upperX, lowerY, lowerZ);
    fReal HHH = getValueAt(upperX, upperY, upperZ);
    fReal HLH = getValueAt(upperX, lowerY, upperZ);

    fReal alphaX = (x - getXCoordAtIndex(lowerX)) * invH;
    fReal alphaY = (y - getYCoordAtIndex(lowerY)) * invH;
    fReal alphaZ = (z - getZCoordAtIndex(lowerZ)) * invH;

    // lower face
    fReal A = Lerp(LLL, LHL, alphaY);
    fReal B = Lerp(LLH, LHH, alphaY);
    fReal C = Lerp(A, B, alphaZ);

    // upper face
    fReal D = Lerp(HLL, HHL, alphaY);
    fReal E = Lerp(HLH, HHH, alphaY);
    fReal F = Lerp(D, E, alphaZ);

    // between faces
    fReal G = Lerp(C, F, alphaX);

    return G;
}

