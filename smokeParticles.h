# pragma once

# include "globalIncludes.h"
# include "smokeQuantity.h"

class SmokeParticles
{
private:
    std::vector<Eigen::Matrix<fReal, 3, 1>> positions;
    fReal h;
public:
    SmokeParticles(fReal h);
    ~SmokeParticles();

    void addParticles(size_t x, size_t y, size_t z, size_t n);
    void updatePositions(SmokeQuantity* u, SmokeQuantity* v, SmokeQuantity* w, fReal deltaT);
    void write_data_bgeo(const std::string& s, const int frame);
};