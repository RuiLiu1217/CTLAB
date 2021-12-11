#ifndef PROJECTION_CONFIGURATION_H
#define PROJECTION_CONFIGURATION_H

class ProjectionConfiguration {
public:
    enum Device {
        SINGLE_CPU = 0,
        MULTI_CPUS = 1,
        SINGLE_GPU = 2,
        MULTI_GPUS = 3
    } device;
    enum Method { // Projection model
        BRANCHLESS = 0,
        VOLUME_RENDERING = 1,
        DOUBLE_PRECISION_BRANCHLESS = 2,
        PSEUDO_DISTANCE_DRIVEN = 3,
        SIDDON = 4,
        BRANCHES_DISTANCE_DRIVEN = 5
    } method;
};

#endif