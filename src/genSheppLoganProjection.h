#ifndef GENSHEPPLOGANPROJECTION_H_
#define GENSHEPPLOGANPROJECTION_H_

#include <vector>
#include "Ellipsoid.h"
std::vector<float> calculateProjection(CTLAB::SheppLoganPhantom& phantom,
	const float sid, const float sdd, const int DNU, std::vector<float> angs, const float col_size, const float col_offset);



#endif
