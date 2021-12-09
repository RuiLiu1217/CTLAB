#ifndef X_RAY_SOURCE_H
#define X_RAY_SOURCE_H
#include <vector>
class XraySource {
private:
    float x0;
    float y0;
    float z0;

    std::vector<float> hangs;
    std::vector<float> hzPos;

public:
    XraySource();
    float getX0() const;
    float getY0() const;
    float getZ0() const;
    float* getAnglePtr();
    float* getZPositions();
    int getPN() const;
};
#endif