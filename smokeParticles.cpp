# include "smokeParticles.h"

SmokeParticles::SmokeParticles(fReal h) : positions(), h(h)
{
}

SmokeParticles::~SmokeParticles()
{
}

void SmokeParticles::addParticles(size_t x, size_t y, size_t z, size_t n){
    for(unsigned int i = 0; i < n; ++i){
        fReal randX = h * static_cast <fReal> (rand()) / static_cast <fReal> (RAND_MAX);
        fReal randY = h * static_cast <fReal> (rand()) / static_cast <fReal> (RAND_MAX);
        fReal randZ = h * static_cast <fReal> (rand()) / static_cast <fReal> (RAND_MAX);

        fReal xVal = x * h + randX;
        fReal yVal = y * h + randY;
        fReal zVal = z * h + randZ;

        Eigen::Matrix<fReal, 3, 1> pos(xVal, yVal, zVal);

        positions.push_back(pos);
    }
}

void SmokeParticles::updatePositions(SmokeQuantity* u, SmokeQuantity* v, SmokeQuantity* w, fReal deltaT)
{
    for(unsigned int i = 0; i < positions.size(); ++i){
        fReal uVel = u->sampleAt(positions[i][0], positions[i][1], positions[i][2]);
        fReal vVel = v->sampleAt(positions[i][0], positions[i][1], positions[i][2]);
        fReal wVel = w->sampleAt(positions[i][0], positions[i][1], positions[i][2]);

        fReal nextX = positions[i][0] + uVel * deltaT;
        fReal nextY = positions[i][1] + vVel * deltaT;
        fReal nextZ = positions[i][2] + wVel * deltaT;

        positions[i][0] = nextX;
        positions[i][1] = nextY;
        positions[i][2] = nextZ;
    }
}

void SmokeParticles::write_data_bgeo(const std::string& s, const int frame)
{
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
            v[k] = 0.0;
        }
    }
    Partio::write(file.c_str(), *parts);
    parts->release();
}