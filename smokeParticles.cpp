# include "smokeParticles.h"

SmokeParticles::SmokeParticles(fReal h) : positions(), h(h)
{
}

SmokeParticles::~SmokeParticles()
{
}

void SmokeParticles::addParticles(size_t x, size_t y, size_t z, size_t n){
    for(unsigned int i = 0; i < n; ++i){
        fReal randTheta = M_2PI * (static_cast <fReal> (rand()) / static_cast <fReal> (RAND_MAX));
        fReal randR = (this->h / 2) * (static_cast <fReal> (rand()) / static_cast <fReal> (RAND_MAX));

        fReal xVal = x * this->h + (this->h / 2) + randR * cos(randTheta);
        fReal yVal = y * this->h + (this->h / 10);
        fReal zVal = z * this->h + (this->h / 2) + randR * sin(randTheta);

        Eigen::Matrix<fReal, 3, 1> pos(xVal, yVal, zVal);
        Eigen::Matrix<fReal, 3, 1> vel(0.0, 0.0, 0.0);

        positions.push_back(pos);
        velocities.push_back(vel);
    }
}

void SmokeParticles::updateParticles(SmokeQuantity* u, SmokeQuantity* v, SmokeQuantity* w, fReal deltaT)
{
    for(unsigned int i = 0; i < positions.size(); ++i){
        fReal uVel = u->sampleAt(positions[i](0, 0), positions[i](1, 0), positions[i](2, 0));
        fReal vVel = v->sampleAt(positions[i](0, 0), positions[i](1, 0), positions[i](2, 0));
        fReal wVel = w->sampleAt(positions[i](0, 0), positions[i](1, 0), positions[i](2, 0));

        velocities[i](0, 0) = uVel;
        velocities[i](1, 0) = vVel;
        velocities[i](2, 0) = wVel;

        fReal nextX = positions[i](0, 0) + uVel * deltaT;
        fReal nextY = positions[i](1, 0) + vVel * deltaT;
        fReal nextZ = positions[i](2, 0) + wVel * deltaT;

        fReal uVelNext = u->sampleAt(nextX, nextY, nextZ);
        fReal vVelNext = v->sampleAt(nextY, nextY, nextZ);
        fReal wVelNext = w->sampleAt(nextX, nextY, nextZ);

        fReal avgUVel = (uVelNext + uVel) / 2.0;
        fReal avgVVel = (vVelNext + vVel) / 2.0;
        fReal avgWVel = (wVelNext + wVel) / 2.0;

        fReal finalX = positions[i](0, 0) + avgUVel * deltaT;
        fReal finalY = positions[i](1, 0) + avgVVel * deltaT;
        fReal finalZ = positions[i](2, 0) + avgWVel * deltaT;

        positions[i](0, 0) = finalX;
        positions[i](1, 0) = finalY;
        positions[i](2, 0) = finalZ;
    }
}

void SmokeParticles::write_data_bgeo(const std::string& s, const int frame)
{
# ifdef PARTIO
    std::string file = s + std::to_string(frame) + ".bgeo";
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    vH = parts->addAttribute("v", Partio::VECTOR, 3);

    for(unsigned int i = 0; i < positions.size(); ++i){
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        float* v = parts->dataWrite<float>(vH, idx);
        for (int k = 0; k < 3; ++k){
            p[k] = positions[i](k, 0);
            v[k] = velocities[i](k, 0);
        }
    }
    Partio::write(file.c_str(), *parts);
    parts->release(); 
# endif
}
