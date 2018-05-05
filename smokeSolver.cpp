# include "smokeSolver.h"

// CONSTRUCTOR / DESTRUCTOR >>>>>>>>>>

SmokeSolver::SmokeSolver(size_t Nx, size_t Ny, size_t Nz, fReal h, fReal alpha, fReal beta, fReal epsilon, fReal sourceSpeed) :
                        Nx(Nx), Ny(Ny), Nz(Nz), h(h), invH(1.0 / h), alpha(alpha), beta(beta), epsilon(epsilon), sourceSpeed(sourceSpeed)
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

    addCenterGridStash("omegaX");
    addCenterGridStash("omegaY");
    addCenterGridStash("omegaZ");

    addCenterGridStash("omegaLength");

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
    force(dt);
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
        for(size_t gridX = 1; gridX < maxX - 1; ++gridX)
        {
            for(size_t gridY = 1; gridY < maxY - 1; ++gridY)
            {
                for(size_t gridZ = 1; gridZ < maxZ - 1; ++gridZ)
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
    buoyancyForce();
    vorticity();
    forceVelocityUpdate(dt);
}

void SmokeSolver::buoyancyForce()
{
    SmokeQuantity* T = attributeTable["T"];
    SmokeQuantity* d = attributeTable["density"];
    GridStash* centerFV = getStashNamed("centerFV");

    // assume an ambient temperature of zero
    // assume an ambient density of zero
    for(size_t i = 1; i < Nx - 1; ++i)
    {
        for(size_t j = 1; j < Ny - 1; ++j)
        {
            for(size_t k = 1; k < Nz - 1; ++k)
            {
                fReal temp = T->getValueAt(i, j, k);
                fReal dens = d->getValueAt(i, j, k);
                fReal buoyancy = -alpha * dens + beta * temp;
                centerFV->writeValueTo(i, j, k, buoyancy);
            }
        }
    }
}

void SmokeSolver::vorticity()
{
    GridStash* uCenter = this->getStashNamed("centerU");
    GridStash* vCenter = this->getStashNamed("centerV");
    GridStash* wCenter = this->getStashNamed("centerW");

    SmokeQuantity* uFace = this->getAttributeNamed("u");
    SmokeQuantity* vFace = this->getAttributeNamed("v");
    SmokeQuantity* wFace = this->getAttributeNamed("w");

    faceToCenterVelocity(uCenter, uFace, 1, 0, 0);
    faceToCenterVelocity(vCenter, vFace, 0, 1, 0);
    faceToCenterVelocity(wCenter, wFace, 0, 0, 1);

    GridStash* omegaX = this->getStashNamed("omegaX");
    GridStash* omegaY = this->getStashNamed("omegaY");
    GridStash* omegaZ = this->getStashNamed("omegaZ");

    fReal centralDiff = 1.0 / (2.0 * h);

    GridStash* omegaLength = this->getStashNamed("omegaLength");

    // calculating and storing angular velocity vector and its length
    for(int i = 1; i < Nx - 1; ++i)
    {
        for(int j = 1; j < Ny - 1; ++j)
        {
            for(int k = 1; k < Nz - 1; ++k)
            {
                fReal x = wCenter->getValueAt(i, j + 1, k) - wCenter->getValueAt(i, j - 1, k);
                x = x - vCenter->getValueAt(i, j, k + 1) + vCenter->getValueAt(i, j, k - 1);
                x = x * centralDiff;
                fReal y = uCenter->getValueAt(i, j, k + 1) - uCenter->getValueAt(i, j, k - 1);
                y = y - wCenter->getValueAt(i + 1, j, k) + wCenter->getValueAt(i - 1, j, k);
                y = y * centralDiff;
                fReal z = vCenter->getValueAt(i + 1, j, k) - vCenter->getValueAt(i - 1, j, k);
                z = z - uCenter->getValueAt(i, j + 1, k) + uCenter->getValueAt(i, j - 1, k);
                z = z * centralDiff;

                omegaX->writeValueTo(i, j, k, x);
                omegaY->writeValueTo(i, j, k, y);
                omegaZ->writeValueTo(i, j, k, z);

                Eigen::Vector3d omega(x, y, z);

                omegaLength->writeValueTo(i, j, k, omega.norm());
            }
        }
    }

    GridStash* centerFU = getStashNamed("centerFU");
    GridStash* centerFV = getStashNamed("centerFV");
    GridStash* centerFW = getStashNamed("centerFW");

    // calculating vorticity force
    for(int i = 1; i < Nx - 1; ++i)
    {
        for(int j = 1; j < Ny - 1; ++j)
        {
            for(int k = 1; k < Nz - 1; ++k)
            {
                fReal ip1 = omegaLength->getValueAt(i + 1, j, k) * centralDiff;
                fReal im1 = omegaLength->getValueAt(i - 1, j, k) * centralDiff;
                fReal jp1 = omegaLength->getValueAt(i, j + 1, k) * centralDiff;
                fReal jm1 = omegaLength->getValueAt(i, j - 1, k) * centralDiff;
                fReal kp1 = omegaLength->getValueAt(i, j, k + 1) * centralDiff;
                fReal km1 = omegaLength->getValueAt(i, j, k - 1) * centralDiff;

                fReal gradOmegaLX = ip1 - im1;
                fReal gradOmegaLY = jp1 - jm1;
                fReal gradOmegaLZ = kp1 - km1;

                Eigen::Vector3d gradOmegaL(gradOmegaLX, gradOmegaLY, gradOmegaLZ);
                fReal gradLength = gradOmegaL.norm();
                gradLength += 1E-20;    // no division by zero

                fReal normalX = gradOmegaLX / gradLength;
                fReal normalY = gradOmegaLY / gradLength;
                fReal normalZ = gradOmegaLZ / gradLength;

                Eigen::Vector3d normal(normalX, normalY, normalZ);

                fReal wX = omegaX->getValueAt(i, j, k);
                fReal wY = omegaY->getValueAt(i, j, k);
                fReal wZ = omegaZ->getValueAt(i, j, k);

                Eigen::Vector3d omega(wX, wY, wZ);

                Eigen::Vector3d vortForce = epsilon * h * normal.cross(omega);

				if (vortForce[2] > 1.0)
				{
					std::cout << vortForce.transpose() << std::endl;
				}

                fReal existingFU = centerFU->getValueAt(i, j, k);
                fReal existingFV = centerFV->getValueAt(i, j, k);
                fReal existingFW = centerFW->getValueAt(i, j, k);

                centerFU->writeValueTo(i, j, k, existingFU + vortForce[0]);
                centerFV->writeValueTo(i, j, k, existingFV + vortForce[1]);
                centerFW->writeValueTo(i, j, k, existingFW + vortForce[2]);
            }
        }
    }
}

void SmokeSolver::forceVelocityUpdate(fReal dt)
{
    GridStash* faceFU = getStashNamed("faceFU");
    GridStash* faceFV = getStashNamed("faceFV");
    GridStash* faceFW = getStashNamed("faceFW");

    GridStash* centerFU = getStashNamed("centerFU");
    GridStash* centerFV = getStashNamed("centerFV");
    GridStash* centerFW = getStashNamed("centerFW");

    centerToFaceStash(faceFU, centerFU, 1, 0, 1);
    centerToFaceStash(faceFV, centerFV, 0, 1, 0);
    centerToFaceStash(faceFW, centerFW, 0, 0, 1);

    SmokeQuantity* u = attributeTable["u"];
    size_t maxX = u->getNx();

    for(size_t i = 1; i < maxX - 1; ++i)
    {
        for(size_t j = 1; j < Ny - 1; ++j)
        {
            for(size_t k = 1; k < Nz - 1; ++k)
            {
                fReal uBeforeUpdate = u->getValueAt(i, j, k);
                fReal deltaU = dt * faceFU->getValueAt(i, j, k);
                u->writeValueTo(i, j, k, uBeforeUpdate + deltaU);
            }
        }
    }

    SmokeQuantity* v = attributeTable["v"];
    size_t maxY = v->getNy();

    for(size_t i = 1; i < Nx - 1; ++i)
    {
        for(size_t j = 1; j < maxY - 1; ++j)
        {
            for(size_t k = 1; k < Nz - 1; ++k)
            {
                fReal vBeforeUpdate = v->getValueAt(i, j, k);
                fReal deltaV = dt * faceFV->getValueAt(i, j, k);
                v->writeValueTo(i, j, k, vBeforeUpdate + deltaV);
            }
        }
    }

    SmokeQuantity* w = attributeTable["w"];
    size_t maxZ = w->getNz();

    for(size_t i = 1; i < Nx - 1; ++i)
    {
        for(size_t j = 1; j < Ny - 1; ++j)
        {
            for(size_t k = 1; k < maxZ - 1; ++k)
            {
                fReal wBeforeUpdate = w->getValueAt(i, j, k);
                fReal deltaW = dt * faceFW->getValueAt(i, j, k);
                w->writeValueTo(i, j, k, wBeforeUpdate + deltaW);
            }
        }
    }

    u->swapBuffer();
    v->swapBuffer();
    w->swapBuffer();
}

void SmokeSolver::projection(fReal dt)
{
    const fReal density = 1.0;
    fReal rhsScale = -this->h * density / dt;
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
    testBVector(b);
    // exit(1);

    // solve for pressure vector
    Eigen::VectorXd pVector(Nx * Ny * Nz);

    Eigen::ConjugateGradient<Eigen::SparseMatrix<fReal>, Eigen::Lower | Eigen::Upper> cg;
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<fReal>, Eigen::Lower, Eigen::IncompleteCholesky<fReal>> cg;
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
    updateVelWithPressure(u, p, scaleP);
    updateVelWithPressure(v, p, scaleP);
    updateVelWithPressure(w, p, scaleP);

    // TEST divergence free condition
    // fillDivergence(b, uSolid, vSolid, wSolid);
    // testDivergenceFree(b);

    u->swapBuffer();
    v->swapBuffer();
    w->swapBuffer();
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
    if(sum > 1E-10){
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
        if(b(i) > 0.5){
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

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

void SmokeSolver::faceToCenterVelocity(GridStash* centerStash, SmokeQuantity *vel, size_t xO, size_t yO, size_t zO)
{
    for(size_t i = 0; i < this->Nx; ++i)
    {
        for (size_t j = 0; j < this->Ny; ++j)
        {
            for (size_t k = 0; k < this->Nz; ++k)
            {
                fReal lowerValue = vel->getValueAt(i, j, k);
                fReal upperValue = vel->getValueAt(i + xO, j + yO, k + zO);
                fReal avg = (lowerValue + upperValue) / 2.0;
                centerStash->writeValueTo(i, j, k, avg);
            }
        }
    }
}

void SmokeSolver::centerToFaceStash(GridStash* faceStash, GridStash* centerStash, size_t xO, size_t yO, size_t zO)
{
    for(size_t i = 0; i < this->Nx - xO; ++i)
    {
        for (size_t j = 0; j < this->Ny - yO; ++j)
        {
            for (size_t k = 0; k < this->Nz - zO; ++k)
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

    // centrally placed source at the bottom of the sim cube
    // helpful to consider a cube with an odd number of subdivisions
    //
    size_t centerX = this->Nx / 2; 
    size_t centerZ = this->Nz / 2;
    size_t bottomY = 2;

    this->setSourceVelocity(v, centerX, bottomY, centerZ, sourceSpeed);
    this->setSourceVelocity(v, centerX, bottomY + 1, centerZ, sourceSpeed);

    // setting temperature and density of sources
    this->setSourceTemperature(centerX, bottomY, centerZ, 1.0);
    this->setSourceDensity(centerX, bottomY, centerZ, 1.0);

    // for a second particle source at the top
    //
    // size_t topY = Ny - 3;
    // this->setSourceVelocity(v, centerX, topY, centerZ, -sourceSpeed);
    // this->setSourceVelocity(v, centerX, topY - 1, centerZ, -sourceSpeed);
    // this->setSourceTemperature(centerX, topY - 1, centerZ, 1.0);
    // this->setSourceDensity(centerX, topY - 1, centerZ, 1.0);

    // for a second and third particle source at the bottom
    //
    // size_t upperLeftX = this->Nx / 4;
    // size_t upperLeftZ = 3 * this->Ny / 4;
    // this->setSourceVelocity(v, upperLeftX, bottomY, upperLeftZ, sourceSpeed);
    // this->setSourceVelocity(v, upperLeftX, bottomY + 1, upperLeftZ, sourceSpeed);
    // this->setSourceTemperature(upperLeftX, bottomY, upperLeftZ, 1.0);
    // this->setSourceDensity(upperLeftX, bottomY, upperLeftZ, 1.0);

    // size_t lowerRightX = 3 * this->Nx / 4;
    // size_t lowerRightZ = this->Ny / 4;
    // this->setSourceVelocity(v, lowerRightX, bottomY, lowerRightZ, sourceSpeed);
    // this->setSourceVelocity(v, lowerRightX, bottomY + 1, lowerRightZ, sourceSpeed);
    // this->setSourceTemperature(lowerRightX, bottomY, lowerRightZ, 1.0);
    // this->setSourceDensity(lowerRightX, bottomY, lowerRightZ, 1.0);

    // for "twin rising flames"
    //
    // size_t uLx = 2 * Nx / 5;
    // size_t uLz = 3 * Nz / 5;
    // size_t uRx = 3 * Nx / 5;
    // size_t uRz = 3 * Nz / 5;
    // size_t lLx = 2 * Nx / 5;
    // size_t lLz = 2 * Nz / 5;
    // size_t lRx = 3 * Nx / 5;
    // size_t lRz = 2 * Nz / 5;
    // size_t bottomY = 2;
    // size_t upperY = Ny / 3;
    // // lower left source
    // this->setSourceVelocity(v, lLx, bottomY, lLz, 2 * sourceSpeed);
    // this->setSourceVelocity(v, lLx, bottomY + 1, lLz, 2 * sourceSpeed);
    // this->setSourceVelocity(u, lLx, bottomY, lLz, sourceSpeed);
    // this->setSourceVelocity(u, lLx + 1, bottomY, lLz, sourceSpeed);
    // this->setSourceVelocity(w, lLx, bottomY, lLz, sourceSpeed / 4.0);
    // this->setSourceDensity(lLx, bottomY, lLz, 1.0);
    // this->setSourceTemperature(lLx, bottomY, lLz, 1.0);
    // // upper right source
    // this->setSourceVelocity(v, uRx, bottomY, uRz, 2 * sourceSpeed);
    // this->setSourceVelocity(v, uRx, bottomY + 1, uRz, 2 * sourceSpeed);
    // this->setSourceVelocity(u, uRx, bottomY, uRz, -sourceSpeed);
    // this->setSourceVelocity(u, uRx, bottomY, uRz + 1, -sourceSpeed);
    // this->setSourceVelocity(w, uRx, bottomY, uRz, -sourceSpeed / 4.0);
    // this->setSourceDensity(uRx, bottomY, uRz, 1.0);
    // this->setSourceTemperature(uRx, bottomY, uRz, 1.0);
}

void SmokeSolver::setSourceVelocity(SmokeQuantity* speed, size_t x, size_t y, size_t z, fReal val)
{    
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

    // GRATE : solver has trouble converging
    // size_t m = 1;
    // size_t n = 1;
    // size_t y = Ny / 2;
    // for(size_t i = 2 * m; i < Nx - 1; ++m)
    // {
    //     for(size_t j = 2 * n; j < Nz - 1; ++n)
    //     {
    //         this->gridTypes[getIndex(i, y, j)] = SOLID;
    //     }
    // }

    // CROSS : works reasonably well
    // size_t left = Nx / 4;
    // size_t right = 3 * Nx / 4;
    // size_t lower = Ny / 3;
    // size_t upper = 2 * Ny / 3;

    // for(size_t j = lower; j <= upper; ++j)
    // {
    //     for(size_t i = left; i < right; ++i)
    //     {
    //         for(size_t k = 0; k < Nz; ++k)
    //         {
    //             this->gridTypes[getIndex(i, j, k)] = SOLID;
    //         }
    //     }
    //     for(size_t k = left; k < right; ++k)
    //     {
    //         for(size_t i = 0; i < Nx; ++i)
    //         {
    //             this->gridTypes[getIndex(i, j, k)] = SOLID;
    //         }
    //     }
    // }
}

// <<<<<<<<<<
// OUTPUT >>>>>>>>>>>>>>>>>>>>>>>>>>>>

void SmokeSolver::write_data_bgeo(const std::string& s, const int frame)
{
# ifdef PARTIO
    std::string file = s + std::to_string(frame) + ".bgeo";
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH, dens;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);
    dens = parts->addAttribute("density", Partio::VECTOR, 1);

    GridStash* uCenter = this->getStashNamed("centerU");
    GridStash* vCenter = this->getStashNamed("centerV");
    GridStash* wCenter = this->getStashNamed("centerW");

    SmokeQuantity* uFace = this->getAttributeNamed("u");
    SmokeQuantity* vFace = this->getAttributeNamed("v");
    SmokeQuantity* wFace = this->getAttributeNamed("w");

    faceToCenterVelocity(uCenter, uFace, 1, 0, 0);
    faceToCenterVelocity(vCenter, vFace, 0, 1, 0);
    faceToCenterVelocity(wCenter, wFace, 0, 0, 1);

    SmokeQuantity* de = this->getAttributeNamed("density");

    for(size_t i = 0; i < this->Nx; ++i)
    {
        for(size_t j = 0; j < this->Ny; ++j)
        {
            for(size_t k = 0; k < this->Nz; ++k)
            {
                int idx = parts->addParticle();
                float* p = parts->dataWrite<float>(posH, idx);
                float* v = parts->dataWrite<float>(vH, idx);
                float* d = parts->dataWrite<float>(dens, idx);

                fReal x, y, z;
                x = (i + 0.5) * h;
                y = (j + 0.5) * h;
                z = (k + 0.5) * h;
                p[0] = x;
                p[1] = y;
                p[2] = z;

                v[0] = uCenter->getValueAt(i, j, k);
                v[1] = vCenter->getValueAt(i, j, k);
                v[2] = wCenter->getValueAt(i, j, k);

                d[0] = de->getValueAt(i, j, k);
            }
        }
    }
    Partio::write(file.c_str(), *parts);
    parts->release();
# endif
}

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

void SmokeSolver::addFaceGridStash(std::string name, size_t xO, size_t yO, size_t zO)
{
    size_t attrNx = this->Nx + xO;
    size_t attrNy = this->Ny + yO;
    size_t attrNz = this->Nz + zO;

    GridStash* ptr = new GridStash(name, attrNx, attrNy, attrNz);
    this->gridStashTable.emplace(std::pair<std::string, GridStash*>(name, ptr));
}

SmokeQuantity* SmokeSolver::getAttributeNamed(std::string name)
{
    return (*this)[name];
}

GridStash* SmokeSolver::getStashNamed(std::string name)
{
    return this->gridStashTable.at(name);
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
    int xIndex = std::floor(x * this->invH);
    int yIndex = std::floor(y * this->invH);
    int zIndex = std::floor(z * this->invH);

    size_t xCell = xIndex < 0 ? 0 : xIndex;
    size_t yCell = yIndex < 0 ? 0 : yIndex;
    size_t zCell = zIndex < 0 ? 0 : zIndex;

    xCell = xCell >= this->Nx ? this->Nx - 1 : xCell;
    yCell = yCell >= this->Ny ? this->Ny - 1 : yCell;
    zCell = zCell >= this->Nz ? this->Nz - 1 : zCell;

    return getGridTypeAt(xCell, yCell, zCell);
}

size_t SmokeSolver::getIndex(size_t x, size_t y, size_t z)
{
    return (z * (this->Nx * this->Ny)) + (y * this->Nx) + x;
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