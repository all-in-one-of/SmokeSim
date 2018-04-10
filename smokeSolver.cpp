# include "smokeSolver.h"

// CONSTRUCTOR / DESTRUCTOR >>>>>>>>>>

SmokeSolver::SmokeSolver(size_t Nx, size_t Ny, size_t Nz, fReal h) :
                        Nx(Nx), Ny(Ny), Nz(Nz), h(h)
{
    addFaceAttr("u", 0.0, 0.5, 0.5);
    addFaceAttr("v", 0.5, 0.0, 0.5);
    addFaceAttr("w", 0.5, 0.5, 0.0);
    addCenterAttr("p", 0.5, 0.5, 0.5);

    this->gridTypes = new gridType[Nx * Ny * Nz];
    memset(reinterpret_cast<void*>(this->gridTypes), SMOKE, Nx * Ny * Nz);
}

SmokeSolver::~SmokeSolver()
{
    for (auto& attr : this->centerAttr)
    {
        delete attr.second;
    }
    for (auto& attr : this->faceAttr)
    {
        delete attr.second;
    }
    delete[] this->gridTypes;
}

// <<<<<<<<<<
// CORE FLUID SOLVER >>>>>>>>>>>>>>>>>

// <<<<<<<<<<
// INITIALIZATION >>>>>>>>>>>>>>>>>>>>

// <<<<<<<<<<
// OUTPUT >>>>>>>>>>>>>>>>>>>>>>>>>>>>

// <<<<<<<<<<
// ACCESS / ATTRIBUTE STORAGE >>>>>>>>

void SmokeSolver::addCenterAttr(std::string name, fReal xOffset, fReal yOffset, fReal zOffset)
{
    size_t attrNx = this->Nx;
    size_t attrNy = this->Ny;
    size_t attrNz = this->Nz;
    
    SmokeQuantity* ptr = new SmokeQuantity(name, attrNx, attrNy, attrNz, this->h, xOffset, yOffset, zOffset);
    this->centerAttr.emplace(std::pair<std::string, SmokeQuantity*>(name, ptr));
}

void SmokeSolver::addFaceAttr(std::string name, fReal xOffset, fReal yOffset, fReal zOffset)
{
    size_t attrNx = this->Nx + 1;
    size_t attrNy = this->Ny + 1;
    size_t attrNz = this->Nz + 1;
    
    SmokeQuantity* ptr = new SmokeQuantity(name, attrNx, attrNy, attrNz, this->h, xOffset, yOffset, zOffset);
    this->faceAttr.emplace(std::pair<std::string, SmokeQuantity*>(name, ptr));
}

SmokeQuantity* SmokeSolver::getAttributeNamed(std::string name)
{
    return (*this)[name];
}

SmokeQuantity* SmokeSolver::operator[](std::string name)
{
    if (centerAttr.find(name) == centerAttr.end())
    {
        return faceAttr.at(name);
    }
    else
    {
        return centerAttr.at(name);
    }
}

// <<<<<<<<<<