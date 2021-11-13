#pragma once


//Calculate Two Arrays inner product
///This function is not used for Users but for calling by system; 
///IT IS A META-PROGRAM
///Author: Rui LIU
///Date: 2013-2-18
///Version: 1.0
template<int Dims, typename T>
class DotProduct
{
public:
	static T result(T* a, T* b)
	{
		return (*a)*(*b) + DotProduct<Dims - 1, T>::result(a + 1, b + 1);
	}
};
template<typename T>
class DotProduct<1, T>
{
public:
	static T result(T* a, T* b)
	{
		return (*a)*(*b);
	}
};

//////////////////////////////////////////////////////////////////////////
/// Call the dot production meta function
template<int Dims, typename T>
inline T innerProd_META(T* a, T* b)
{
	return DotProduct<Dims, T>::result(a, b);
}



//Add two vector together
//This function is not used for Users but for calling by system;
//It is a META-program
//Author: Rui LIU
//Date: 2013-8-9
//Version: 1.0
template<int Dims, typename T>
class AddVector
{
public:
	static void result(T* a, T* b, T* c)
	{
		(*c) = (*a) + (*b);
		AddVector<Dims - 1, T>::result(a + 1, b + 1);
	}
};
template<typename T>
class AddVector<1, T>
{
public:
	static T result(T* a, T* b, T* c)
	{
		(*c) = (*a) + (*b);
	}
};

template<int Dims, typename T>
inline void add_META(T* a, T* b, T* c)
{
	return AddVector<Dims, T>::result(a, b, c);
}


//Sub vector a from b together
//This function is not used for Users but for calling by system;
//It is a META-program
//Author: Rui LIU
//Date: 2013-8-9
//Version: 1.0
template<int Dims, typename T>
class SubVector
{
public:
	static void result(T* a, T* b, T* c)
	{
		(*c) = (*a) - (*b);
		SubVector<Dims - 1, T>::result(a + 1, b + 1);
	}
};
template<typename T>
class SubVector<1, T>
{
public:
	static T result(T* a, T* b, T* c)
	{
		(*c) = (*a) - (*b);
	}
};

template<int Dims, typename T>
inline void sub_META(T* a, T* b, T* c)
{
	return SubVector<Dims, T>::result(a, b, c);
}



///This metaprogram is used to calculate the p power for each element and then add them together
//It is a META-PROGRAM
//Author: Rui Liu
//Date: 2013-08-13
//Version: 1.0
template<int Dims, typename T>
class powerVector
{
public:
	static T result(T* vec, T p)
	{
		return std::pow((*vec), p) + powerVector<Dims - 1, T>::result(vec + 1, p);
	}
};
template<typename T>
class powerVector<1, T>
{
public:
	static T result(T *vec, T p)
	{
		return std::pow((*vec), p);
	}
};



///This meta function solute the problem that a matrix multiply with
/// a vector
///Author: Rui Liu
/// Use the meta programming method to generate the 
/// cosine and sine series given a series of theta
template<int Dims, typename T>
class CosSinTheta
{
public:
	static void result(T* cosTheta, T* sinTheta, T* theta)
	{
		(*cosTheta) = cos((*theta));
		(*sinTheta) = sin((*theta));
		CosSinTheta<Dims - 1, T>::result(cosTheta + 1, sinTheta + 1, theta + 1);
	}
};
template<typename T>
class CosSinTheta<1, T>
{
public:
	static T result(T* cosTheta, T* sinTheta, T* theta)
	{
		(*cosTheta) = cos((*theta));
		(*sinTheta) = sin((*theta));
	}
};
