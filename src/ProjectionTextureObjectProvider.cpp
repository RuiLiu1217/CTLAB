#include "ProjectionTextureObjectProvider.h"

ProjectionTextureObjectProvider::~ProjectionTextureObjectProvider() {
    destory();
}

cudaTextureObject_t ProjectionTextureObjectProvider::getTexObj() const {
    return texObj;
}
cudaArray* ProjectionTextureObjectProvider::getArray() const {
    return d_prjArray;
}

void ProjectionTextureObjectProvider::destory() {
    cudaDestroyTextureObject(texObj);
    cudaFreeArray(d_array);
}