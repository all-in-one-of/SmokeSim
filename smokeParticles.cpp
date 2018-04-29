# include "smokeParticles.h"

SmokeParticles::SmokeParticles(fReal h) : positions(), h(h)
{
}

SmokeParticles::~SmokeParticles()
{
}

void SmokeParticles::addParticles(size_t x, size_t y, size_t z, size_t n){
    for(unsigned int i = 0; i < n; ++i){
        fReal randTheta = M_2PI * static_cast <fReal> (rand()) / static_cast <fReal> (RAND_MAX);
        fReal randR = (h / 2) * static_cast <fReal> (rand()) / static_cast <fReal> (RAND_MAX);

        fReal xVal = x * h + (h / 2) + randR * cos(randTheta);
        fReal yVal = y * h;
        fReal zVal = z * h + (h / 2) + randR * sin(randTheta);

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

        positions[i](0, 0) = nextX;
        positions[i](1, 0) = nextY;
        positions[i](2, 0) = nextZ;
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
            v[k] = velocities[i](k, 0);
        }
    }
    Partio::write(file.c_str(), *parts);
    parts->release();
}
