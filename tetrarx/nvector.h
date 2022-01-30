#pragma once

#include<iostream>
#define TINY 0.0000001f

using namespace std ;

class nVector
{
public:
	int n;
	float *x;      //координаты
	void *tag;

	nVector(int n_)
	{
		n = n_;
		x = new float[n];
		for(int i = 0; i < n; ++i) x[i] = 0.0f;
	};
		
	nVector(int n_, float val)
	{
		n = n_;
		x = new float[n];
		for(int i = 0; i < n; ++i) x[i] = val;
	};

	nVector()
	{
		n = 0;
		x = nullptr;
		tag = nullptr;
	};
	
	void initnew(nVector v)
	{
		if(n > 0)
			delete [] x;
		n = v.n;
		x = new float[n];
		for(int i = 0; i < n; ++i) x[i] = v.x[i];
	}

	void initnew(int n_)
	{
		if(n > 0)
			delete [] x;
		n = n_;
		x = new float[n];
		for(int i = 0; i < n; ++i) x[i] = 0.0f;
	};

	inline nVector operator+ (const nVector& b) const
	{
		nVector a(n);
		for(int i = 0; i < n; ++i) a.x[i] = x[i] + b.x[i];
		return a;
	};

	inline nVector operator+ (const float b) const
	{
		nVector a(n);
		for(int i = 0; i < n; ++i) a.x[i] = x[i] + b;
		return a;
	};

	inline nVector operator- () const
	{
		nVector a(n);
		for(int i = 0; i < n; ++i) a.x[i] = -x[i];
		return a;
	};

	inline nVector operator- (const nVector& b) const
	{
		nVector a(n);
		for(int i = 0; i < n; ++i) a.x[i] = x[i] - b.x[i];
		return a;
	};

	inline nVector operator* (const float b) const
	{
		nVector a(n);
		for(int i = 0; i < n; ++i) a.x[i] = x[i] * b;
		return a;
	};

	inline const float operator* (const nVector & b) const                 //скал€рное пр-е
	{
		float sum = 0.0f;
		for(int i = 0; i < n; ++i) sum += x[i] * b.x[i];
		return sum;
	};

	inline nVector & operator= (const nVector &b)
	{
		if(this == &b || b.n == 0)
			return *this;

		if(n == 0) {
			n = b.n;
			x = new float[n];
		}
		else if(n != b.n)
		{
			delete [] x;
			n = b.n;
			x = new float[n];
		}

		for(int i = 0; i < n; ++i) x[i] = b.x[i];
		return *this;
	};

	inline nVector & operator= (const float b)
	{
		for(int i = 0; i < n; ++i) x[i] = b;
		return *this;
	};

	inline float dr(nVector b)
	{
		float sum = 0.0f;
		float dif;
		
		for(int i = 0; i < n; ++i)
		{
			dif = b.x[i] - x[i];
			sum += dif * dif;
		}
		return sqrtf(sum);
	};

	inline float sqr(void) const
	{
		float sum = 0.0f;
		for(int i = 0; i < n; ++i) sum += x[i] * x[i];
		return sum;
	};

	inline float norm(void) const
	{
        return sqrtf(sqr());
    };

	inline void stabilize_length() 
	{
		float sum = 0.0f;
		for(int i = 0; i < n; ++i) sum += fabs(x[i]);
     
	   if(sum > 0.0f){ for(int i = 0; i < n; ++i) x[i] /= sum; }
	   else {          for(int i = 0; i < n; ++i) x[i] = 0.0f; }
	};

		
	inline void scale(const float s) const
	{
		for(int i = 0; i < n; ++i) x[i] = x[i] * s;
	};

	inline void normalize(void) {
		float m = norm();
        if( m < TINY )
        {
            stabilize_length();
			m = norm();
        }

		if(m > 0.0f)
			scale( 1.0f/m );
    };


	//~nVector()
	//{
	//	if(n > 0 && x != nullptr)
	//		delete [] x;
	//}


};

	
inline const nVector operator+ (const float c, const nVector & b)
{
	nVector a(b.n);
	for(int i = 0; i < b.n; ++i) a.x[i] = c + b.x[i];
	return a;
};

inline const nVector operator- (const float c, const nVector & b)
{
	nVector a(b.n);
	for(int i = 0; i < b.n; ++i) a.x[i] = c - b.x[i];
	return a;
};

inline const nVector operator* (const float c, const nVector & b)
{
	nVector a(b.n);
	for(int i = 0; i < b.n; ++i) a.x[i] = c * b.x[i];
	return a;
};

