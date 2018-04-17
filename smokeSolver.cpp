# include "smokeSolver.h"

// CONSTRUCTOR / DESTRUCTOR >>>>>>>>>>

SmokeSolver::SmokeSolver(size_t Nx, size_t Ny, size_t Nz, fReal h) :
                        Nx(Nx), Ny(Ny), Nz(Nz), h(h)
{
    addFaceAttr("u", -0.5, 0.0, 0.0);
    addFaceAttr("v", 0.0, -0.5, 0.0);
    addFaceAttr("w", 0.0, 0.0, -0.5);
    addCenterAttr("density", 0.0, 0.0, 0.0);
    addCenterAttr("T", 0.0, 0.0, 0.0);

    this->gridTypes = new gridType[Nx * Ny * Nz];
    memset(reinterpret_cast<void*>(this->gridTypes), SMOKE, Nx * Ny * Nz);
}

SmokeSolver::~SmokeSolver()
{
    for (auto& attr : this->attributeTable)
    {
        delete attr.second;
    }
    for (auto& attr : this->gridStashTable)
    {
        delete attr.second;
    }
    delete[] this->gridTypes;
}

// <<<<<<<<<<
// CORE FLUID SOLVER >>>>>>>>>>>>>>>>>

void SmokeSolver::step(fReal dt)
{
    // advection(dt);
    // bodyForce(dt);
    // projection(dt);
}

void KaminoSolver::advection(fReal dt)
{
    for (auto quantity : this->attributeTable)
    {
        SmokeQuantity* attr = quantity.second;
        size_t maxX = attr->getNx();
        size_t maxY = attr->getNy();
        size_t maxZ = attr->getNz();
        for(size_t gridX = 0; gridX < maxX; ++gridX)
        {
            for(size_t gridY = 0; gridY < maxY; ++gridY)
            {
                for(size_t gridZ = 0; gridZ < maxZ; ++gridZ)
                {
                    fReal gX = attr->getXCoordAtIndex(gridX);
                    fReal gY = attr->getYCoordAtIndex(gridY);
                    fReal gZ = attr->getZCoordAtIndex(gridZ);

                    fReal uG = (*this)["u"]->sampleAt(gX, gY, gZ);
                    fReal vG = (*this)["v"]->sampleAt(gX, gY, gZ);
                    fReal wG = (*this)["w"]->sampleAt(gX, gY, gZ);

                    fReal midX = gX - 0.5 * dt * uG;
                    fReal midY = gY - 0.5 * dt * vG;
                    fReal midZ = gZ - 0.5 * dt * wG;

                    fReal uMid = (*this)["u"]->sampleAt(midX, midY, midZ);
                    fReal vMid = (*this)["v"]->sampleAt(midX, midY, midZ);
                    fReal wMid = (*this)["w"]->sampleAt(midX, midY, midZ);

                    fReal pX = gX - dt * uMid;
                    fReal pY = gY - dt * vMid;
                    fReal pZ = gZ - dt * wMid;
                    
                    fReal advectedVal = attr->sampleAt(pX, pY, pZ);
                    attr->writeValueTo(gridX, gridY, gridZ, advectedVal);
                }
            }
        }
    }
    this->swapBuffers();
}

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
    //this->centerAttr.emplace(std::pair<std::string, SmokeQuantity*>(name, ptr));
    this->attributeTable.emplace(std::pair<std::string, SmokeQuantity*>(name, ptr));
}

void SmokeSolver::addFaceAttr(std::string name, fReal xOffset, fReal yOffset, fReal zOffset)
{
    size_t attrNx = this->Nx + 1;
    size_t attrNy = this->Ny + 1;
    size_t attrNz = this->Nz + 1;
    
    SmokeQuantity* ptr = new SmokeQuantity(name, attrNx, attrNy, attrNz, this->h, xOffset, yOffset, zOffset);
    //this->faceAttr.emplace(std::pair<std::string, SmokeQuantity*>(name, ptr));
    this->attributeTable.emplace(std::pair<std::string, SmokeQuantity*>(name, ptr));
}

void SmokeSolver::addGridStash(std::string name)
{
    size_t attrNx = this->Nx;
    size_t attrNy = this->Ny;
    size_t attrNz = this->Nz;

    gridStash* ptr = new GridStash(name, attrNx, attrNy, attrNz);
    this->gridStashTable.emplace(std::pair<std::string, GridStash*>(name, ptr));
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

gridType SmokeSolver::getGridTypeAt(size_t x, size_t y, size_t z)
{
    return gridTypes[getIndex(x, y, z)];
}

void SmokeSolver::swapBuffers()
{
    for (auto quantity : this->attributeTable)
    {
        quantity.second->swapBuffer();
    }
}

void SmokeSolver::swapFaceBuffers()
{
    for (auto quantity : this->faceAttr)
    {
        quantity.second->swapBuffer();
    }
}

void KaminoSolver::swapCenterBuffers()
{
    for (auto quantity : this->centerAttr)
    {
        quantity.second->swapBuffer();
    }
}

// <<<<<<<<<<