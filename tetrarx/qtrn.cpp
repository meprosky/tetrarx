#define TINY 0.000001f
#include<math.h>

class Quaternion
{
public:
	float x, y, z, w;
	
	Quaternion(float x = 0.0f, float y = 0.0f, float z = 0.0f, float w = 0.0f);
	
	inline void ident(void)
	 {
		 x = 0.0f; y = 0.0f; z = 0.0f; w = 1.0f; 
	 }

	inline void scale(float s)
	{
		x *= s; y *= s; z *= s; 
	}
	
	inline float norm(void) const
	{
		return x*x + y*y + z*z + w*w;
	}
	
	inline float magnitude(void) const
	{
		return sqrtf( norm() );
	}

	inline Quaternion inverted(void) const
	{
        float in = 1.0f / norm();
        return Quaternion( -x*in, -y*in, -z*in, w*in );
    } 

	inline void invert(void)
	 {
		 float in = 1.0f / norm();
		 x *= -in; y *= -in; z *= -in; w *= in;
	 }

	/**
      stabilize quaternion length within 1 - 1/4
      this operation is a lot faster than normalization
      and preserve length goes to 0 or infinity
    */
    inline void stabilize_length()
	{
		float cs = fabs(x) + fabs(y) + fabs(z) + fabs(w);
		if( cs > 0.0f )
		{
			x /= cs; y /= cs; z /= cs; w /= cs; 
		}
		else
		{
			ident();
		}
    }

	 /**
        scale quaternion that its norm goes to 1 
        the appearing of 0 magnitude or near is a error,
        so we can be sure that can divide by magnitude
    */

    inline void normalize(void)
	{
        float m = magnitude();
        if( m < TINY )
        {
            stabilize_length();
            m = magnitude();
        }
        scale( 1.0f/m );
    }


	inline const Quaternion& operator+=(const Quaternion& q)
	{
        x += q.x; y += q.y; z += q.z; w += q.w;
		return *this;
	}

	inline const Quaternion& operator-=(const Quaternion& q)
	{
        x -= q.x; y -= q.y; z -= q.z; w -= q.w;
		return *this;
	}
	 

	/*
	  Regular mult
      q( cross(v,v') + wv' + w'v, ww' - dot(v,v') )
    */
	/*
    inline  Quaternion operator*(const  Quaternion& q)const {
        return  Quaternion( w*q.x + x*q.w + y*q.z - z*q.y,
                            w*q.y + y*q.w + z*q.x - x*q.z,
                            w*q.z + z*q.w + x*q.y - y*q.x,
                            w*q.w - x*q.x - y*q.y - z*q.z);
        // 16 multiplications    12 addidtions    0 variables
    }
	*/

	
	inline  Quaternion operator*( const Quaternion& q) const
	{
        float t0 = (x-y)*(q.y-q.x);
        float t1 = (w+z)*(q.w+q.z);
        float t2 = (w-z)*(q.y+q.x);
        float t3 = (x+y)*(q.w-q.z);
        float t4 = (x-z)*(q.z-q.y);
        float t5 = (x+z)*(q.z+q.y);
        float t6 = (w+y)*(q.w-q.x);
        float t7 = (w-y)*(q.w+q.x);

        float t8 = t5 + t6 + t7;
        float t9 = (t4 + t8)*0.5f;
        return Quaternion ( t3+t9-t6,
                            t2+t9-t7,
                            t1+t9-t8,
                            t0+t9-t5 );

        // 9 multiplications    27  addidtions    8 variables
        // but of couse we can clean 4 variables
		
	}
	


	
	inline const Quaternion& operator*=(const Quaternion& q)
	{
		Quaternion res = *this * q;
        x = res.x;
		y = res.y;
		z = res.z;
		w = res.w;
		
		//set((*this)*q);  // have no optimization here
    }
	

	



/*
int Quaternion::normalize()
{
	float a2 = angle_rad / 2;
	float d = sinf(a2);

	

	xn = x * d;
	yn = y * d;
	zn = z * d;
	wn = cosf(a2);

	d = sqrtf(xn*xn + yn*yn + zn*zn + wn*wn);
	
	if(d == 0.0f)
		return -1;

	d = 1.0f / d;

	xn = xn * d;
	yn = yn * d;
	zn = zn * d;
	wn = wn * d;

	return 1;
}
*/
	



}; //end of Class





