#include "XraySource.h"
XraySource::XraySource() {
    x0 = 0.0f;
    y0 = 500.0f;
    z0 = 0.0f;

    int PN = 720;
    hangs.resize(PN);
    hzPos.resize(PN);

    for (int i = 0; i < PN; ++i) {
        hangs[i] = 2.0 * 3.141592653589793 / PN * i;
        hzPos[i] = 0.0;
    }
}
float XraySource::getX0() const {
    return x0;
}

float XraySource::getY0() const {
    return y0;
}

float XraySource::getZ0() const {
    return z0;
}

float* XraySource::getAnglePtr() {
    return &(hangs[0]);
}

float* XraySource::getZPositions() {
    return &(hzPos[0]);
}

int XraySource::getPN() const {
    return hangs.size();
}