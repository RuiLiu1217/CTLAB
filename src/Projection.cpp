#include "Projection.h"
Projection::Projection() {
    data.resize(888 * 512 * 720);
}

int Projection::getDNU() const {
    return detector.getDNU();
}
int Projection::getDNV() const {
    return detector.getDNV();
}

int Projection::getPN() const {
    return xraySource.getPN();
}

float* Projection::getDataPtr() {
    return &(data[0]);
}

float Projection::getX0() const {
    return xraySource.getX0();
}
float Projection::getY0() const {
    return xraySource.getY0();
}

float Projection::getZ0() const {
    return xraySource.getZ0();
}

float* Projection::getXdsPtr() {
    return detector.getXdsPtr();
}

float* Projection::getYdsPtr() {
    return detector.getYdsPtr();
}

float* Projection::getZdsPtr() {
    return detector.getZdsPtr();
}

float* Projection::getAnglePtr() {
    return xraySource.getAnglePtr();
}

float* Projection::getZPosPtr() {
    return xraySource.getZPositions();
}