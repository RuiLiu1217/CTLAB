/*!
 * \file Ellipsoid.h
 * \brief Definition of Ellipsoid for 3D Shepp-Logan phantom
 * Wake Forest University Health Sciences
 * \author Rui Liu
 * \date Jun 14, 2016
 * \version 1.0
 */
#ifndef ELLIPSOID_H_
#define ELLIPSOID_H_

#include <cmath>
#include <iostream>
#include <vector>
/// \brief CTLAB is the main namespace for our library
namespace CTLAB
{
	/// \brief Basic datatype with two single floating values
	struct float2
	{
		float x;
		float y;

		float2() :x(0), y(0){}
		float2(const float a, const float b) :x(a), y(b){}
		float2(const CTLAB::float2& rgh) :x(rgh.x), y(rgh.y){}

		float2& operator=(const float2& rgh)
		{
			if (this != &rgh)
			{
				x = rgh.x;
				y = rgh.y;
			}
			return *this;
		}

		float2 operator+=(const float2& rhs)
		{
			x += rhs.x;
			y += rhs.y;
			return *this;
		}

		friend float2 operator+(float2 lhs, const float2& rhs)
		{
			lhs += rhs;
			return lhs;
		}

		float2 operator-=(const float2& rhs)
		{
			x -= rhs.x;
			y -= rhs.y;
			return *this;
		}

		friend float2 operator-(float2 lhs, const float2& rhs)
		{
			lhs -= rhs;
			return lhs;
		}

		/// \brief rotate the point with a given angle arc
		/// \return new point position
		float2 rot(const float arc)
		{
			float2 res;
			float cosT = cosf(arc);
			float sinT = sinf(arc);
			res.x = x * cosT - y * sinT;
			res.y = x * sinT + y * cosT;
			return res;
		}

		float2 invRot(const float arc)
		{
			return rot(-arc);
		}

		float length()
		{
			return sqrtf(x*x + y*y);
		}

		float2 normalize()
		{
			float2 res;
			float len = length();
			res.x = x / len;
			res.y = y / len;
			return res;
		}


	};
	
	/// \brief Ellipsoid
	class Ellipsoid
	{
	public:
		Ellipsoid(){};
		virtual ~Ellipsoid(){};
	private:
		float2 center;
		float2 axis;
		float rotAng; // rotation angle with respect to x-axis in arc
		float intensity; // intensity of the  ellipsoid.
	public:
		Ellipsoid(float2 cnt, float2 axs, float rot, float intens) :center(cnt), axis(axs), rotAng(rot), intensity(intens){}

		Ellipsoid(const Ellipsoid& rgh) :center(rgh.getCenter()), axis(rgh.getAxis()), rotAng(rgh.getRotAng()), intensity(rgh.getIntensity()){}

		Ellipsoid& operator=(const Ellipsoid& rgh)
		{
			if (this != &rgh)
			{
				center = rgh.getCenter();
				axis = rgh.getAxis();
				rotAng = rgh.getRotAng();
				intensity = rgh.getIntensity();
			}
			return *this;
		}


		const float2 getCenter() const {
			return center;
		}
		bool setCenter(const float2& v)
		{
			center = v;
			return true;
		}

		bool setAxis(const float2& x)
		{
			axis = x;
			return true;
		}
		bool setCenter(const float x, const float y)
		{
			center.x = x;
			center.y = y;
			return true;
		}

		bool setAxis(const float x, const float y)
		{
			axis.x = x;
			axis.y = y;
			return true;
		}

		const float2 getAxis() const{
			return axis;
		}
		const float getCenterX() const
		{
			return center.x;
		}
		const float getCenterY() const
		{
			return center.y;
		}
		const float getAxisX() const
		{
			return axis.x;
		}
		const float getAxisY() const
		{
			return axis.y;
		}

		const float getRotAng() const
		{
			return rotAng;
		}

		bool setRotAng(const float ang)
		{
			rotAng = ang;
			return true;
		}

		const float getIntensity() const
		{
			return intensity;
		}

		bool setIntensity(const float intens)
		{
			intensity = intens;
			return true;
		}

		float intersectionLength(float2 start, float2 end)
		{
			//Move the line and ellipsoid that start and end point both intersect with the centered ellipsoids.
			start -= center;
			end -= center;



			//rotate the line
			float2 nStart = start.invRot(rotAng);
			float2 nEnd = end.invRot(rotAng);

			float2 dir = nEnd - nStart;
			dir = dir.normalize();


			const float a = getAxisX();
			const float b = getAxisY();
			const float invA2 = 1.0 / (a * a);
			const float invB2 = 1.0 / (b * b);

			const float A = (dir.x * dir.x * invA2 + dir.y * dir.y * invB2);
			const float B = 2.0 * (dir.x * nStart.x * invA2 + dir.y * nStart.y * invB2);
			const float C = nStart.x * nStart.x * invA2 + nStart.y * nStart.y * invB2 - 1.0;

			const float Delta = B * B - 4.0 * A * C;
			if (Delta <= 0)
			{
				return 0;
			}
			else
			{
				const float t1 = (-B + sqrtf(Delta)) / (2.0 * A);
				const float t2 = (-B - sqrtf(Delta)) / (2.0 * A);

				return fabs(t1 - t2);
			}
			return 0;
		}

	};

	/// \brief 3D Shepp Logan phantom
	class SheppLoganPhantom
	{

	public:
		std::vector<Ellipsoid> ele;
		SheppLoganPhantom()
		{
			ele.resize(10);
			//		[  1   .69   .92    0     0     0]
			ele[0].setIntensity(1.0);
			ele[0].setAxis(0.69 / 2.0, 0.92 / 2.0);
			ele[0].setCenter(0.0, 0.0);
			ele[0].setRotAng(0);
			//		 -.8  .6624 .8740   0  -.0184   0
			ele[1].setIntensity(-0.8);
			ele[1].setAxis(0.6624 / 2.0, 0.8740 / 2.0);
			ele[1].setCenter(0.0, -0.0184 / 2);
			ele[1].setRotAng(0);
			//		 -.2  .1100 .3100  .22    0    -18
			ele[2].setIntensity(-0.2);
			ele[2].setAxis(0.11 / 2.0, 0.31 / 2.0);
			ele[2].setCenter(0.22 / 2, 0.0);
			ele[2].setRotAng(-18.0 / 180.0 * 3.141592653589793);
			//		 -.2  .1600 .4100 -.22    0     18
			ele[3].setIntensity(-0.2);
			ele[3].setAxis(0.16 / 2.0, 0.41 / 2.0);
			ele[3].setCenter(-0.22 / 2, 0.0);
			ele[3].setRotAng(18.0 / 180.0 * 3.141592653589793);
			//		  .1  .2100 .2500   0    .35    0
			ele[4].setIntensity(0.1);
			ele[4].setAxis(0.21 / 2.0, 0.25 / 2.0);
			ele[4].setCenter(0, 0.35 / 2);
			ele[4].setRotAng(0);
			//		  .1  .0460 .0460   0    .1     0
			ele[5].setIntensity(0.1);
			ele[5].setAxis(0.046 / 2.0, 0.046 / 2.0);
			ele[5].setCenter(0, 0.1 / 2);
			ele[5].setRotAng(0);

			//		  .1  .0460 .0460   0   -.1     0
			ele[6].setIntensity(0.1);
			ele[6].setAxis(0.046 / 2.0, 0.046 / 2.0);
			ele[6].setCenter(0, -0.1 / 2);
			ele[6].setRotAng(0);
			//		  .1  .0460 .0230 -.08  -.605   0
			ele[7].setIntensity(0.1);
			ele[7].setAxis(0.046 / 2.0, 0.023 / 2.0);
			ele[7].setCenter(-0.08 / 2, -0.605 / 2);
			ele[7].setRotAng(0);
			//		  .1  .0230 .0230   0   -.606   0
			ele[8].setIntensity(0.1);
			ele[8].setAxis(0.023 / 2.0, 0.023 / 2.0);
			ele[8].setCenter(0, -0.606 / 2);
			ele[8].setRotAng(0);
			//		  .1  .0230 .0460  .06  -.605   0   ];
			ele[9].setIntensity(0.1);
			ele[9].setAxis(0.023 / 2.0, 0.046 / 2.0);
			ele[9].setCenter(0.06 / 2, -0.605 / 2);
			ele[9].setRotAng(0);
		}

		~SheppLoganPhantom()
		{
			ele.clear();
		}
	};


} /* namespace CTLAB */
#endif
