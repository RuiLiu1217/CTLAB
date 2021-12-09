#ifndef PROJECTION_H_
#define PROJECTION_H_
#include <vector>
#include "Detector.h"
#include "XraySource.h"
class Projection {
private:
    std::vector<float> data;

    XraySource xraySource;
    Detector detector;
public:
    Projection();
    int getDNU() const;
    int getDNV() const;
    int getPN() const;
    float* getDataPtr();
    
    float getX0() const;
    float getY0() const;
    float getZ0() const;
    float* getXdsPtr();
    float* getYdsPtr();
    float* getZdsPtr();
    float* getAnglePtr();
    float* getZPosPtr();

};

#endif