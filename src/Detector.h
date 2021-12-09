#ifndef DETECTOR_H
#define DETECTOR_H
#include<vector>
class Detector {
private:
    std::vector<float> xds;
    std::vector<float> yds;
    std::vector<float> zds;
public:
    Detector();
    float* getXdsPtr();
    float* getYdsPtr();
    float* getZdsPtr();
    int getDNU() const;
    int getDNV() const;
};
#endif