/*! 
 * \file genSheppLoganProjection.h
 * \brief Generate the 2D Fan beam analytical projection for SVD paper
 * \author Rui Liu
 * \version 1.0
 * \date Jun. 1, 2014
 */
#ifndef GENSHEPPLOGANPROJECTION_H_
#define GENSHEPPLOGANPROJECTION_H_

#include <vector>
#include "Ellipsoid.h"
/// \brief calculate the analytical projection 
/// \param phantom The object of 3D phantom, which can be modified but the modified function has not been implemented.
/// \param sid source to iso-center distance
/// \param sdd source to detector distance
/// \param DNU number of detector cell along channel direction
/// \param angs all views
/// \param col_size detector size 
/// \param col_offset offset of the detector along column direction (channel direction)
std::vector<float> calculateProjection(CTLAB::SheppLoganPhantom& phantom,
	const float sid, const float sdd, const int DNU, std::vector<float> angs, const float col_size, const float col_offset);

#endif
