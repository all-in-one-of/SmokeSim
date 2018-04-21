# include "smokeSolver.h"

// CONSTRUCTOR / DESTRUCTOR >>>>>>>>>>

SmokeSolver::SmokeSolver(size_t Nx, size_t Ny, size_t Nz, fReal h) :
                        Nx(Nx), Ny(Ny), Nz(Nz), h(h), invH(1.0 / h)
{
    addFaceAttr("u", -0.5, 0.0, 0.0);
    addFaceAttr("v", 0.0, -0.5, 0.0);
    addFaceAttr("w", 0.0, 0.0, -0.5);

    addCenterAttr("density", 0.0, 0.0, 0.0);
    addCenterAttr("T", 0.0, 0.0, 0.0);

    addGridStash("pressure", 0);

    addGridStash("avgU", 0);
    addGridStash("avgV", 0);
    addGridStash("avgW", 0);

    this->gridTypes = new gridType[Nx * Ny * Nz];
    memset(reinterpret_cast<void*>(this->gridTypes), SMOKE, Nx * Ny * Nz);
    initializeBoundary();
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
    advection(dt);
    // bodyForce(dt);
    // projection(dt);
}

void SmokeSolver::advection(fReal dt)
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

                    gridType g = getGridTypeAt(pX, pY, pZ);

                    if(g == SMOKE){
                        fReal advectedVal = attr->sampleAt(pX, pY, pZ);
                        attr->writeValueTo(gridX, gridY, gridZ, advectedVal);
                    }
                    else{   // point to advect from is inside SOLID grid cell
                        fReal advectedVal = 0.0;
                        attr->writeValueTo(gridX, gridY, gridZ, advectedVal);
                    }
                }
            }
        }
    }
    this->swapBuffers();
}

void SmokeSolver::averageAttribute(GridStash* avgStash, SmokeQuantity* attr)
{
    size_t maxX = attr->getNx() - 1;
    size_t maxY = attr->getNy() - 1;
    size_t maxZ = attr->getNz() - 1;

    for(size_t i = 0; i < maxX; ++i)
    {
        for (size_t j = 0; j < maxY; ++j)
        {
            for (size_t k = 0; k < maxZ; ++k)
            {
                fReal lowerValue = attr->getValueAt(i, j, k);
                fReal upperValue = attr->getValueAt(i + 1, j + 1, k + 1);
                fReal avg = (lowerValue + upperValue) / 2.0;
                avgStash->writeValueTo(i, j, k, avg);
            }
        }
    }
}

// <<<<<<<<<<
// INITIALIZATION >>>>>>>>>>>>>>>>>>>>

void SmokeSolver::setSourceVelocity(size_t x, size_t y, size_t z, fReal val){
    SmokeQuantity* v = (*this)["v"];
    v->setValueAt(x, y, z, val);
    v->setValueAt(x, y + 1, z, val);
}

void SmokeSolver::initializeBoundary(){
    // vw planes
    for(size_t j = 0; j < this->Ny; ++j){
        for(size_t k = 0; k < this->Nz; ++k){
            // left face
            this->gridTypes[getIndex(0, j, k)] = SOLID;
            // right face
            this->gridTypes[getIndex(Nx - 1, j, k)] = SOLID;
        }
    }
    // uv planes
    for(size_t i = 0; i < this->Nx; ++i){
        for(size_t j = 0; j < this->Ny; ++j){
            // back face
            this->gridTypes[getIndex(i, j, 0)] = SOLID;
            // front face
            this->gridTypes[getIndex(i, j, Nz - 1)] = SOLID;
        }
    }
    // uw planes
    for(size_t i = 0; i < this->Nx; ++i){
        for(size_t k = 0; k < this->Nz; ++k){
            // top face
            this->gridTypes[getIndex(i, 0, k)] = SOLID;
            // bottom face
            this->gridTypes[getIndex(i, Ny - 1, k)] = SOLID;
        }
    }
}

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
    size_t attrNx = this->Nx;
    size_t attrNy = this->Ny;
    size_t attrNz = this->Nz;

    if(name == "u"){
        attrNx++;
    }
    if(name == "v"){
        attrNy++;
    }
    if(name == "w"){
        attrNz++;
    }
    
    SmokeQuantity* ptr = new SmokeQuantity(name, attrNx, attrNy, attrNz, this->h, xOffset, yOffset, zOffset);
    //this->faceAttr.emplace(std::pair<std::string, SmokeQuantity*>(name, ptr));
    this->attributeTable.emplace(std::pair<std::string, SmokeQuantity*>(name, ptr));
}

void SmokeSolver::addGridStash(std::string name, int offset)
{
    size_t attrNx = this->Nx + offset;
    size_t attrNy = this->Ny + offset;
    size_t attrNz = this->Nz + offset;

    GridStash* ptr = new GridStash(name, attrNx, attrNy, attrNz);
    this->gridStashTable.emplace(std::pair<std::string, GridStash*>(name, ptr));
}

SmokeQuantity* SmokeSolver::getAttributeNamed(std::string name)
{
    return (*this)[name];
}

SmokeQuantity* SmokeSolver::operator[](std::string name)
{
    return attributeTable.at(name);
    // if (centerAttr.find(name) == centerAttr.end())
    // {
    //     return faceAttr.at(name);
    // }
    // else
    // {
    //     return centerAttr.at(name);
    // }
}

gridType SmokeSolver::getGridTypeAt(size_t x, size_t y, size_t z)
{
    return gridTypes[getIndex(x, y, z)];
}

gridType SmokeSolver::getGridTypeAt(fReal x, fReal y, fReal z)
{
    // const fReal epsilon = 1E-10;
    int xIndex = std::floor((x + 0.5) * invH);
    int yIndex = std::floor((y + 0.5) * invH);
    int zIndex = std::floor((z + 0.5) * invH);

    size_t xCell = xIndex < 0 ? 0 : xIndex;
    size_t yCell = yIndex < 0 ? 0 : yIndex;
    size_t zCell = zIndex < 0 ? 0 : zIndex;

    xCell = xCell >= Nx ? Nx - 1 : xCell;
    yCell = yCell >= Ny ? Ny - 1 : yCell;
    zCell = zCell >= Nz ? Nz - 1 : zCell;

    return getGridTypeAt(xCell, yCell, zCell);
}

size_t SmokeSolver::getIndex(size_t x, size_t y, size_t z)
{
    return z * (Nx * Ny) + y * (Nx) + x;
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

void SmokeSolver::swapCenterBuffers()
{
    for (auto quantity : this->centerAttr)
    {
        quantity.second->swapBuffer();
    }
}

// <<<<<<<<<<