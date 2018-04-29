# include "smokeSolver.h"

// CONSTRUCTOR / DESTRUCTOR >>>>>>>>>>

SmokeSolver::SmokeSolver(size_t Nx, size_t Ny, size_t Nz, fReal h) :
                        Nx(Nx), Ny(Ny), Nz(Nz), h(h), invH(1.0 / h)
{
    addFaceAttr("u", 0.0, 0.5, 0.5);
    addFaceAttr("v", 0.5, 0.0, 0.5);
    addFaceAttr("w", 0.5, 0.5, 0.0);

    addCenterAttr("density", 0.5, 0.5, 0.5);
    addCenterAttr("T", 0.5, 0.5, 0.5);

    addCenterGridStash("pressure");
    addCenterGridStash("buoyancy");

    addCenterGridStash("centerU");
    addCenterGridStash("centerV");
    addCenterGridStash("centerW");

    addCenterGridStash("centerFU");
    addCenterGridStash("centerFV");
    addCenterGridStash("centerFW");
    
    addFaceGridStash("faceFU", 1, 0, 0);
    addFaceGridStash("faceFV", 0, 1, 0);
    addFaceGridStash("faceFW", 0, 0, 1);

    this->gridTypes = new gridType[Nx * Ny * Nz];
    memset(reinterpret_cast<void*>(this->gridTypes), SMOKE, Nx * Ny * Nz);
    initializeBoundary();

    precomputeLaplacian();
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
    setSources();
    advection(dt);
    // force(dt);
    projection(dt);
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

void SmokeSolver::force(fReal dt)
{
    fReal alpha = 0.08;
    fReal beta = 1.0;

    SmokeQuantity* T = attributeTable["T"];
    SmokeQuantity* d = attributeTable["density"];
    GridStash* centerFV = getStashNamed("centerFV");

    // buoyancy force
    for(size_t i = 0; i < Nx; ++i)
    {
        for(size_t j = 0; j < Ny; ++j)
        {
            for(size_t k = 0; k < Nz; ++k)
            {
                fReal temp = T->getValueAt(i, j, k);
                fReal dens = d->getValueAt(i, j, k);
                fReal buoyancy = -alpha * dens + beta * temp;
                centerFV->writeValueTo(i, j, k, buoyancy);
            }
        }
    }

    GridStash* faceFV = getStashNamed("faceFV");
    centerToFaceStash(faceFV, centerFV, 0, 1, 0);

    SmokeQuantity* v = attributeTable["v"];
    size_t maxY = v->getNy();

    for(size_t i = 0; i < Nx; ++i)
    {
        for(size_t j = 0; j < maxY; ++j)
        {
            for(size_t k = 0; k < Nz; ++k)
            {
                fReal vBeforeUpdate = v->getValueAt(i, j, k);
                fReal deltaV = dt * faceFV->getValueAt(i, j, k);
                v->writeValueTo(i, j, k, vBeforeUpdate + deltaV + 500.0);
            }
        }
    }

    v->swapBuffer();
}

void SmokeSolver::projection(fReal dt)
{
    const fReal density = 1.0;
    fReal rhsScale = -h * density / dt;
    fReal scaleP = 1.0 / rhsScale;

    Eigen::VectorXd b(Nx * Ny * Nz);
    b.setZero();

    SmokeQuantity* u = attributeTable["u"];
    SmokeQuantity* v = attributeTable["v"];
    SmokeQuantity* w = attributeTable["w"];

    fReal uSolid = 0.0;
    fReal vSolid = 0.0;
    fReal wSolid = 0.0;

    // filling rhs divergence values
    fillDivergence(b, uSolid, vSolid, wSolid);

    b = b * rhsScale;

    // TEST that b vector sums to zero
    // testBVector(b);
    // exit(1);

    // solve for pressure vector
    Eigen::VectorXd pVector(Nx * Ny * Nz);

    Eigen::ConjugateGradient<Eigen::SparseMatrix<fReal>, Eigen::Lower | Eigen::Upper> cg;
    cg.compute(Laplacian);
    pVector = cg.solve(b);

    GridStash* p = getStashNamed("pressure");

    // populate updated pressure values
    for(size_t i = 1; i < Nx - 1; ++i)
    {
        for(size_t j = 1; j < Ny - 1; ++j)
        {
            for(size_t k = 1; k < Nz - 1; ++k)
            {
                p->writeValueTo(i, j, k, pVector(getIndex(i, j, k)));
            }
        }
    }

    // update velocities
    for (auto quantity : this->faceAttr){
        updateVelWithPressure(quantity.second, p, scaleP);
    }

    u->swapBuffer();
    v->swapBuffer();
    w->swapBuffer();

    // TEST divergence free condition
    // fillDivergence(b, uSolid, vSolid, wSolid);
    // testDivergenceFree(b);
}

void SmokeSolver::updateVelWithPressure(SmokeQuantity* speed, GridStash* p, fReal scaleP)
{
    size_t xMax = speed->getNx();
    size_t yMax = speed->getNy();
    size_t zMax = speed->getNz();

    size_t uOffSet = 0;
    size_t vOffSet = 0;
    size_t wOffSet = 0;

    if(speed->getName() == "u"){
        uOffSet = 1;
    }
    if(speed->getName() == "v"){
        vOffSet = 1;
    }
    if(speed->getName() == "w"){
        wOffSet = 1;
    }

    for(size_t i = 1; i < xMax - 1; ++i)
    {
        for(size_t j = 1; j < yMax - 1; ++j)
        {
            for(size_t k = 1; k < zMax - 1; ++k)
            {
                fReal speedBeforeUpdate = speed->getValueAt(i, j, k);

                if(getGridTypeAt(i, j, k) == SMOKE && getGridTypeAt(i - uOffSet, j - vOffSet, k - wOffSet) == SMOKE){
                    fReal deltaPressure = p->getValueAt(i, j, k) - p->getValueAt(i - uOffSet, j - vOffSet, k - wOffSet);
                    fReal delta = scaleP * deltaPressure;
                    speed->writeValueTo(i, j, k, speedBeforeUpdate + delta);
                }
                else{
                    speed->writeValueTo(i, j, k, 0.0);
                }
            }
        }
    }
}

void SmokeSolver::testBVector(Eigen::VectorXd& b)
{
    fReal sum = 0.0;
    size_t size = Nx * Ny * Nz;
    for(int i = 0; i < size; ++i){
        sum += b(i);
    }
    std::cout << sum << std::endl;
    if(sum > 0.0){
        std::cout << "error: b vector does not sum to zero" << std::endl;
    }
}

void SmokeSolver::fillDivergence(Eigen::VectorXd& b, fReal uSolid, fReal vSolid, fReal wSolid)
{
    SmokeQuantity* u = attributeTable["u"];
    SmokeQuantity* v = attributeTable["v"];
    SmokeQuantity* w = attributeTable["w"];

    for(size_t i = 1; i < Nx - 1; ++i)
    {
        for(size_t j = 1; j < Ny - 1; ++j)
        {
            for(size_t k = 1; k < Nz - 1; ++k)
            {
                if(getGridTypeAt(i, j, k) == SOLID){
                    continue;
                }

                fReal uPlus, uMinus, vPlus, vMinus, wPlus, wMinus;
                size_t ipoot = i + 1;
                size_t imoot = i;
                size_t jpoot = j + 1;
                size_t jmoot = j;
                size_t kpoot = k + 1;
                size_t kmoot = k;

                if(getGridTypeAt(i + 1, j, k) == SMOKE){
                    uPlus = u->getValueAt(ipoot, j, k);
                }
                else{
                    uPlus = uSolid;
                }
                if(getGridTypeAt(i - 1, j, k) == SMOKE){
                    uMinus = u->getValueAt(imoot, j, k);
                }
                else{
                    uMinus = uSolid;
                }
                if(getGridTypeAt(i, j + 1, k) == SMOKE){
                    vPlus = v->getValueAt(i, jpoot, k);
                }
                else{
                    vPlus = vSolid;
                }
                if(getGridTypeAt(i, j - 1, k) == SMOKE){
                    vMinus = v->getValueAt(i, jmoot, k);
                }
                else{
                    vMinus = vSolid;
                }
                if(getGridTypeAt(i, j, k + 1) == SMOKE){
                    wPlus = w->getValueAt(i, j, kpoot);
                }
                else{
                    wPlus = wSolid;
                }
                if(getGridTypeAt(i, j, k - 1) == SMOKE){
                    wMinus = w->getValueAt(i, j, kmoot);
                }
                else{
                    wMinus = wSolid;
                }
                b(getIndex(i, j, k)) = (uPlus - uMinus) + (vPlus - vMinus) + (wPlus - wMinus);
            }
        }
    }
}

void SmokeSolver::testDivergenceFree(Eigen::VectorXd& b)
{
    size_t size = Nx * Ny * Nz;
    for(size_t i = 0; i < size; ++i)
    {
        if(b(i) > 1E-10){
            std::cout << "error: divergence of " << b(i) << std::endl;
        }
    }
}

void SmokeSolver::precomputeLaplacian()
{
    Laplacian = Eigen::SparseMatrix<fReal>(Nx * Ny * Nz, Nx * Ny * Nz);
    Laplacian.setZero();

    // interior of boundary
    for(size_t i = 1; i < Nx - 1; ++i)
    {
        for(size_t j = 1; j < Ny - 1; ++j)
        {
            for(size_t k = 1; k < Nz - 1; ++k)
            {
                size_t numXNeighbors = 0;
                size_t numYNeighbors = 0;
                size_t numZNeighbors = 0;

                size_t ip1 = (i + 1);
                size_t im1 = (i - 1);
                size_t jp1 = (j + 1);
                size_t jm1 = (j - 1);
                size_t kp1 = (k + 1);
                size_t km1 = (k - 1);

                size_t rowNumber = getIndex(i, j, k);

                if(getGridTypeAt(i, j, k) == SOLID){
                   continue; 
                }

                if(getGridTypeAt(ip1, j, k) == SMOKE){
                    Laplacian.coeffRef(rowNumber, getIndex(ip1, j, k)) = -1;
                    numXNeighbors++;
                }
                if(getGridTypeAt(im1, j, k) == SMOKE){
                    Laplacian.coeffRef(rowNumber, getIndex(im1, j, k)) = -1;
                    numXNeighbors++;
                }
                if(getGridTypeAt(i, jp1, k) == SMOKE){
                    Laplacian.coeffRef(rowNumber, getIndex(i, jp1, k)) = -1;
                    numYNeighbors++;
                }
                if(getGridTypeAt(i, jm1, k) == SMOKE){
                    Laplacian.coeffRef(rowNumber, getIndex(i, jm1, k)) = -1;
                    numYNeighbors++;
                }
                if(getGridTypeAt(i, j, kp1) == SMOKE){
                    Laplacian.coeffRef(rowNumber, getIndex(i, j, kp1)) = -1;
                    numZNeighbors++;
                }
                if(getGridTypeAt(i, j, km1) == SMOKE){
                    Laplacian.coeffRef(rowNumber, getIndex(i, j, km1)) = -1;
                    numZNeighbors++;
                }
                Laplacian.coeffRef(rowNumber, getIndex(i, j, k)) = numXNeighbors + numYNeighbors + numZNeighbors;
            }
        }
    }

    // TEST that Laplacian is symmetric and rows sum to zero
    // testLaplacian();
    // exit(1);
}

void SmokeSolver::testLaplacian()
{
    size_t size = Nx * Ny * Nz;
    Eigen::VectorXd test(size);
    for(size_t i = 0; i < size; ++i){
        test(i) = 1.0;
    }  
    Eigen::VectorXd zero(size);
    zero = Laplacian * test;
    for(size_t i = 0; i < size; ++i){
        if(zero(i) > 0){
            std::cout << "error: Laplacian row doesn't sum to zero" << std::endl;
        }
    }
    Eigen::SparseMatrix<fReal> LaplacianT(size, size);
    LaplacianT = Laplacian.transpose();
    zero = LaplacianT * test;
    for(size_t i = 0; i < size; ++i){
        if(zero(i) > 0){
            std::cout << "error: Laplacian isn't symmetric" << std::endl;
        }
    }
}

void SmokeSolver::faceToCenterVelocity(GridStash* centerStash, SmokeQuantity* vel, int xO, int yO, int zO)
{
    size_t maxX = vel->getNx() - xO;
    size_t maxY = vel->getNy() - yO;
    size_t maxZ = vel->getNz() - zO;

    for(size_t i = 0; i < maxX; ++i)
    {
        for (size_t j = 0; j < maxY; ++j)
        {
            for (size_t k = 0; k < maxZ; ++k)
            {
                fReal lowerValue = vel->getValueAt(i, j, k);
                fReal upperValue = vel->getValueAt(i + xO, j + yO, k + zO);
                fReal avg = (lowerValue + upperValue) / 2.0;
                centerStash->writeValueTo(i, j, k, avg);
            }
        }
    }
}

void SmokeSolver::centerToFaceStash(GridStash* faceStash, GridStash* centerStash, int xO, int yO, int zO)
{
    size_t maxX = centerStash->getNx() - xO;
    size_t maxY = centerStash->getNy() - yO;
    size_t maxZ = centerStash->getNz() - zO;

    for(size_t i = 0; i < maxX; ++i)
    {
        for (size_t j = 0; j < maxY; ++j)
        {
            for (size_t k = 0; k < maxZ; ++k)
            {
                fReal lowerValue = centerStash->getValueAt(i, j, k);
                fReal upperValue = centerStash->getValueAt(i + xO, j + yO, k + zO);
                fReal avg = (lowerValue + upperValue) / 2.0;
                faceStash->writeValueTo(i + xO, j + yO, k + zO, avg);
            }
        }
    }
}

// <<<<<<<<<<
// INITIALIZATION >>>>>>>>>>>>>>>>>>>>

void SmokeSolver::setSources()
{
    SmokeQuantity* u = (*this)["u"];
    SmokeQuantity* v = (*this)["v"];
    SmokeQuantity* w = (*this)["w"];
    // if we assume a centrally placed source at the bottom of the sim cube
    // helpful to consider a cube with an odd number of subdivisions
    size_t centerX = Nx / 2; 
    size_t centerZ = Nz / 2;
    size_t bottomY = 1;
    this->setSourceVelocity(v, centerX, bottomY, centerZ, 2.0);
    this->setSourceVelocity(v, centerX, bottomY + 1, centerZ, 2.0);

    // setting temperature and density of sources
    this->setSourceTemperature(centerX, bottomY, centerZ, 1.0);
    this->setSourceTemperature(centerX, bottomY + 1, centerZ, 1.0);
    this->setSourceDensity(centerX, bottomY, centerZ, 1.0);
    this->setSourceDensity(centerX, bottomY + 1, centerZ, 1.0);
}

void SmokeSolver::setSourceVelocity(SmokeQuantity* speed, size_t x, size_t y, size_t z, fReal val){
    
    speed->setValueAt(x, y, z, val);
}

void SmokeSolver::setSourceTemperature(size_t x, size_t y, size_t z, fReal val){
    SmokeQuantity* T = (*this)["T"];
    T->setValueAt(x, y, z, val);
}

void SmokeSolver::setSourceDensity(size_t x, size_t y, size_t z, fReal val){
    SmokeQuantity* d = (*this)["density"];
    d->setValueAt(x, y, z, val);
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
    this->attributeTable.emplace(std::pair<std::string, SmokeQuantity*>(name, ptr));
}

void SmokeSolver::addFaceAttr(std::string name, fReal xOffset, fReal yOffset, fReal zOffset)
{
    size_t attrNx = this->Nx;
    size_t attrNy = this->Ny;
    size_t attrNz = this->Nz;

    if(name == "u"){
        attrNx = attrNx + 1;
    }
    if(name == "v"){
        attrNy = attrNy + 1;
    }
    if(name == "w"){
        attrNz = attrNz + 1;
    }
    
    SmokeQuantity* ptr = new SmokeQuantity(name, attrNx, attrNy, attrNz, this->h, xOffset, yOffset, zOffset);
    this->attributeTable.emplace(std::pair<std::string, SmokeQuantity*>(name, ptr));
}

void SmokeSolver::addCenterGridStash(std::string name)
{
    size_t attrNx = this->Nx;
    size_t attrNy = this->Ny;
    size_t attrNz = this->Nz;

    GridStash* ptr = new GridStash(name, attrNx, attrNy, attrNz);
    this->gridStashTable.emplace(std::pair<std::string, GridStash*>(name, ptr));
}

void SmokeSolver::addFaceGridStash(std::string name, int xO, int yO, int zO)
{
    size_t attrNx = this->Nx + 1;
    size_t attrNy = this->Ny + 1;
    size_t attrNz = this->Nz + 1;

    GridStash* ptr = new GridStash(name, attrNx, attrNy, attrNz);
    this->gridStashTable.emplace(std::pair<std::string, GridStash*>(name, ptr));
}

SmokeQuantity* SmokeSolver::getAttributeNamed(std::string name)
{
    return (*this)[name];
}

GridStash* SmokeSolver::getStashNamed(std::string name)
{
    return gridStashTable.at(name);
}

SmokeQuantity* SmokeSolver::operator[](std::string name)
{
    return attributeTable.at(name);
}

gridType SmokeSolver::getGridTypeAt(size_t x, size_t y, size_t z)
{
    return gridTypes[getIndex(x, y, z)];
}

gridType SmokeSolver::getGridTypeAt(fReal x, fReal y, fReal z)
{
    int xIndex = std::floor(x * invH);
    int yIndex = std::floor(y * invH);
    int zIndex = std::floor(z * invH);

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