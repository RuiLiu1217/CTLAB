#include "Detector.h"

Detector::Detector() {
    int DNU = 888;
    int DNV = 512;

    float y0 = 500.0f;
    xds.resize(DNU);
    yds.resize(DNU);
    zds.resize(DNV);
    
    float R = 300.0;
    float theta = 3.141592653589793 / 4.5;

    float delta_theta = theta / DNU;

    for (int i = 0; i < DNU; ++i) {
        xds[i] = sin(-theta / 2.0 + delta_theta * i) * R;
        yds[i] = -cos(-theta / 2.0 + delta_theta * i) * R + y0;
    }

    for (int i = 0; i < DNV; ++i) {
        zds[i] = ( - 256.0f + i * 1.0f) / 2.4f;
    }
}
float* Detector::getXdsPtr() {
    return &(xds[0]);
}

float* Detector::getYdsPtr() {
    return &(yds[0]);
}

float* Detector::getZdsPtr() {
    return &(zds[0]);
}

int Detector::getDNU() const {
    return xds.size();
}

int Detector::getDNV() const {
    return zds.size();
}