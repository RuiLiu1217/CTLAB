#ifndef VOLUME_TO_TEXTURES_H
#define VOLUME_TO_TEXTURES_H


class VolumeToTextures {
private:
    cudaArray* d_volumeArray1;
    cudaArray* d_volumeArray2;
    cudaTextureObject_t texObj1;
    cudaTextureObject_t texObj2;

public:
    
    // Assume that vol is in host memory
    VolumeToTextures(float* vol, int XN, int YN, int ZN, cudaTextureAddressMode textureAddressMode);
    ~VolumeToTextures();

    cudaTextureObject_t getTextureZXY() const;
    cudaTextureObject_t getTextureZYX() const;
};

class VolumeToOneTexture {
private:
    cudaArray* d_volumeArray;
    cudaTextureObject_t texObj;

public:
    VolumeToOneTexture(float* vol, int XN, int YN, int ZN, cudaTextureAddressMode textureAddressMode);
    ~VolumeToOneTexture();
    cudaTextureObject_t getTexture() const;
};

#endif