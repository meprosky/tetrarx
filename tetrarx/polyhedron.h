#pragma once

#include<iostream>
#include<math.h>
#include<limits>
#include"param.h"
#include"nvector.h"


#pragma comment(lib, "rdrand64.lib") // апп. ген. сл. чисел

using namespace std ;

int _rdrand32_step(unsigned int *);


inline unsigned int rand32()  //random 0..RAND_MAX32 (4294967295)
{
	unsigned int x;
    _rdrand32_step(&x);
	return x;
}
inline unsigned int rand32ui(const unsigned int min, const unsigned int max) //[min..max)   т.е. сл.ч.  < max,  min > 0
{
	unsigned int x = rand32();
	
	if (RAND_MAX32 == x) return rand32ui(min, max);

	unsigned int 
		range       = max - min,
		remainder   = RAND_MAX32 % range,
		bucket      = RAND_MAX32 / range;
	
	if (x > RAND_MAX32 - remainder) 
		return rand32ui(min, max);
	else
		return min + x / bucket;
}
inline float rand32fi(const float min, const float max) 
{
	return (float)rand32() / RAND_MAX32 * (max - min) + min;
}

//enum Mendeleev{H=1, C=6, N=7, O=8, Na=11, Mg=12, Al=13, Si=14, K=19, Ca=20, Cr=24, Mn=25, Fe=26, Cu=29, Zn=30, Mo=42, W=74};
//const int Mendeleev_n = 17;

class MendeleevPT
{
public:
	int n;
	int atomnum;
	string name;

	MendeleevPT(int n_, int a_, string name_) : n(n_), atomnum(a_), name(name_) {};

};

enum PT {unknown = 0, Si=1, O=2};

static MendeleevPT periodictable[] = { MendeleevPT(0,0,""), MendeleevPT(1,14,"Si"), MendeleevPT(2,8,"O")};
static const int num_sorts = 3;
static float atomscollen[num_sorts][num_sorts];




class Vector
{
public:
	float x, y, z, w;      //координаты

	Vector( float x_ = 0.0f, float y_ = 0.0f, float z_ = 0.0f, float w_ = 0.0f) : x(x_), y(y_), z(z_), w(w_) {};
		
	inline Vector operator+ (const Vector& b) const
	{
		return Vector(x + b.x, y + b.y, z + b.z);
	}

	inline Vector operator+ (const float b) const
	{
		return Vector(x+b, y+b, z+b);
	}

	inline Vector operator- () const
	{
		return Vector(-x, -y, -z);
	}

	inline Vector operator- (const Vector& b) const
	{
		return Vector(x - b.x, y - b.y, z - b.z);
	}

	inline Vector operator* (const float b) const
	{
		return Vector(x*b, y*b, z*b);
	}

	inline const Vector operator* (const Vector & b) const                 //вект. произв.
	{
	return Vector(
		y*b.z - z*b.y, 
		z*b.x - x*b.z, 
		x*b.y - y*b.x
		);
	};

	inline const float operator& (const Vector & b) const                //скаляр. произв.
	{
		return  x*b.x +  y*b.y + z*b.z;
	};

	inline Vector & operator= (const Vector &b)
	{
		if(this == &b)
			return *this;
		x = b.x; y = b.y; z = b.z; w = b.w;
		return *this;
	};

	inline Vector & operator= (const float b)
	{
		x = b; y = b; z = b; w = b;
		return *this;
	};

	inline void scale3(float s)
	{
		x *= s; y *= s; z *= s; 
	}
	
	inline void scale4(float s)
	{
		x *= s; y *= s; z *= s; w *= s; 
	}

	inline float dr(Vector b)
	{
		float x = b.x - x;
		float y = b.y - y;
		float z = b.z - z;

		return sqrtf(x*x + y*y + z*z);
	};

	inline float norm3(void) const
	{
		return x*x + y*y + z*z;
	};

	inline float norm4(void) const
	{
		return x*x + y*y + z*z + w*w;
	};

	inline float magnitude3(void) const
	{
        return sqrtf(norm3());
    };

	inline float magnitude4(void) const
	{
		return sqrtf(norm4());
    };

	inline void stabilize_length3() 
	{
       float cs = (float)(fabs(x) + fabs(y) + fabs(z));
	   if( cs > 0.0f ){ x/=cs;  y/=cs;  z/=cs;  }
	   else           { x=0.0f; y=0.0f; z=0.0f; }
	}

	inline void stabilize_length4() 
	{
       float cs = (float)(fabs(x) + fabs(y) + fabs(z) + fabs(w));
	   if( cs > 0.0f ){ x/=cs;  y/=cs;  z/=cs;  w/=cs;}
	   else           { x=0.0f; y=0.0f; z=0.0f; w=0.0f;}
	}

	inline void normalize3(void) {
        float m = magnitude3();
        if( m < TINY )
        {
            stabilize_length3();
            m = magnitude3();
        }
		if(m > 0.0f)
			scale3( 1.0f/m );
    }

	inline void normalize4(void) {
        float m = magnitude4();
        if( m < TINY )
        {
            stabilize_length4();
            m = magnitude4();
        }
        scale4( 1.0f/m );
    }

	inline void rotate(Vector q, float norm)  //быстрое вращение
	{
		
		float wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;
        x2 = q.x + q.x; y2 = q.y + q.y; z2 = q.z + q.z;

        xx = q.x * x2;   xy = q.x * y2;   xz = q.x * z2;
        yy = q.y * y2;   yz = q.y * z2;   zz = q.z * z2;
        wx = q.w * x2;   wy = q.w * y2;   wz = q.w * z2;

        x2 = x  - x * (yy + zz) + y * (xy - wz) + z * (xz + wy);
        y2 = y  + x * (xy + wz) - y * (xx + zz) + z * (yz - wx);
		z2 = z  + x * (xz - wy) + y * (yz + wx) - z * (xx + yy);

		x = x2 / norm; y = y2 / norm; z = z2 / norm;
	}

	inline void rotate_my(Vector q, float norm)
	{
		
		Vector v(x,y,z);
		
		float aa, bb, cc, ww, ab, cw, ac, bw, bc, aw;
		float aabb,  abab, acac, cwcw, bcbc, awaw, bwbw, aa_bb, cc_ww;

		aa = q.x*q.x; bb = q.y*q.y; cc = q.z*q.z; ww = q.w*q.w;
		ab = q.x*q.y; ac = q.x*q.z; aw = q.x*q.w; 
		bc = q.y*q.z; bw = q.y*q.w;	cw = q.z*q.w;	

		aabb = aa + bb; abab = ab + ab; acac = ac + ac; cwcw = cw + cw; bcbc = bc + bc;
		awaw = aw + aw; bwbw = bw + bw; aa_bb = aa - bb; cc_ww = cc - ww;;
		
		x = ( aa_bb-cc_ww)*v.x + (abab-cwcw)*v.y + (acac+bwbw)*v.z;
		y = (-aa_bb-cc_ww)*v.y + (abab+cwcw)*v.x + (bcbc-awaw)*v.z;
		z = (-aa-bb+cc+ww)*v.z + (acac-bwbw)*v.x + (bcbc+awaw)*v.y;

		x /= norm; y /= norm; z /= norm;

	}

	inline void rotate_slow(Vector q, float norm)
	{

		Vector v(x,y,z);
		
		Vector      qq(v.x * q.w + v.z * q.y - v.y * q.z,
                       v.y * q.w + v.x * q.z - v.z * q.x,
                       v.z * q.w + v.y * q.x - v.x * q.y,
                       v.x * q.x + v.y * q.y + v.z * q.z);

        x = q.w * qq.x + q.x * qq.w + q.y * qq.z - q.z * qq.y;
        y = q.w * qq.y + q.y * qq.w + q.z * qq.x - q.x * qq.z;
        z = q.w * qq.z + q.z * qq.w + q.x * qq.y - q.y * qq.x;
		
		x /= norm; y /= norm; z /= norm;

	
	}

};

inline const Vector operator+ (const float a, const Vector & b)
{
	return Vector(a+b.x, a+b.y, a+b.z);
};

inline const Vector operator- (const float a, const Vector & b)
{
	return Vector(a-b.x, a-b.y, a-b.z);
};

inline const Vector operator* (const float a, const Vector & b)
{
	return Vector(a*b.x, a*b.y, a*b.z);
};

inline float dr(const Vector & a, const Vector & b)
{
	return (b - a).magnitude3();
}

inline float anglebetweentwovectors(const Vector & a, const Vector & b)
{
	Vector va = a;
    Vector vb = b;
	va.normalize3();
	vb.normalize3();
    
	float sp = va & vb;

	if(fabsf(sp) > 1.0f) sp = (sp < 0) ? -1.0f : 1.0f;
	
	float ang = acosf(sp);
    
	return ang;

	/*


	float ab = a.magnitude3() * b.magnitude3();
	
	if(ab < TINY)
		return 0.0f;

	float ac = (a & b) / ab;
	
	if(ac > 1.0f)
		ac = 1.0f;

	ac = acosf(ac);

	return ac; //acosf( (a & b) / ab);
	*/
}

class Atom
{
public:	
	Vector lc, gc;      //локальная и глобальная координата атома
	
	int anum, sort;     //атомный номер, сорт атома
	string name;        //название в период.табл. Менделеева

	int ida;            //глобальный идентификатор

	int deleted;
	int connected;
	int allowtoconnect;
	
	int nbonds;
	Atom **bonds;

	Atom & operator = (const Atom & other)
    {
        if (this != &other) // защита от неправильного самоприсваивания
        {
            // 1: выделяем "новую" память и копируем элементы
            Atom **new_bonds;
			new_bonds = new Atom* [other.nbonds];

			for(int i = 0; i < other.nbonds; ++i)
				new_bonds[i] = other.bonds[i];

			// 2: освобождаем "старую" память
			delete [] bonds;
 
            // 3: присваиваем значения в "новой" памяти объекту
            bonds  = new_bonds;
			
			nbonds = other.nbonds;

			lc = other.lc;
			gc = other.gc;
	
			anum = other.anum;
			sort = other.sort;
			name = other.name;
			
			ida            = other.ida;

			deleted        = other.deleted;
			connected      = other.connected;
			allowtoconnect = other.allowtoconnect;
        }
        
		// по соглашению всегда возвращаем *this
        return *this;
    }

	Atom(float x_ = 0.0f, float y_ = 0.0f, float z_ = 0.0f, int element_ = 0, int ida_ = -1) :
		lc(x_, y_, z_), gc(x_, y_, z_), ida(ida_), nbonds(0), bonds(0), connected(0), allowtoconnect(1), deleted(0)  
	{
		name = periodictable[element_].name;
		anum = periodictable[element_].atomnum;
		sort = periodictable[element_].n;
	};

	inline void addbond(Atom *a)
	{
		Atom **b;
		b = new Atom* [nbonds+1];

		for(int i = 0; i < nbonds; ++i)
		{
			b[i] = bonds[i];
		}

		b[nbonds] = a;

		if(nbonds > 0)
			delete [] bonds;
		
		++nbonds;
		bonds = b;
	}

	inline void addbond2(Atom *a)
	{
		addbond(a);
		(*a).addbond(this);
	}

	inline void deleteinbondedatoms()
	{
		for(int i = 0; i < nbonds; ++i)
		{
			(*bonds[i]).deletebond(this);
		}

		deleteallbonds();
	}

	inline void deleteallbonds()
	{
		delete [] bonds;
		nbonds = 0;
		bonds = 0;
	}

	inline void replacebonds(Atom *a)
	{
		for(int i = 0; i < nbonds; ++i)
		{
			int id = (*bonds[i]).findbond(this);
			if(id >= 0)
				(*bonds[i]).bonds[id] = a;
		}

		deleteallbonds();

	}

	inline int deletebond(Atom *a)
	{
		Atom **b;
		int j = 0;
		b = new Atom* [nbonds-1];

		for(int i = 0; i < nbonds; ++i)
		{
			if(bonds[i] == a) continue;

			if(j + 1 == nbonds - 1)
			{
				delete [] b;
				return -1;
			}
			b[j++] = bonds[i];
		}
		delete [] bonds;
		--nbonds;
		bonds = b;
		return 0;
	}
	
	inline int findbond(Atom *a)
	{
		for(int i = 0; i < nbonds; ++i)
		{
			if(bonds[i] == a)
				return i;
		}
		return -1;
	}

};


class Polyhedron
{
public:
	static int npall;

	int num_sort;

	int state;			//состояние
	
	int state2;			//состояние

	int idp;            //глобальный идентификатор
	int nvertex;        
	int na;             //число атомов в a

	float L;			//длина грани правильного многогранника
	float rfactor;      //R-фактор для отбора конфигураций
	
	Polyhedron *topoly;
	int fromvertex, tovertex;

	nVector point;      //вектор состояний

	Vector xe, ye, ze;	//единичные вектора осей x,y,z внутри многогр.

	Vector r0;			//относительно начала координат

	Vector qtrn;		//квартернион вращения
	Vector qtrn_r0;		//r0 квартерниона
	float  qtrn_norm;

	Atom *a;   //атомы

	Polyhedron & operator=(const Polyhedron & other)
    {
		if (this != &other) // защита от неправильного самоприсваивания
        {
            // 1: выделяем "новую" память и копируем элементы
			Atom *new_atoms = new Atom[other.na];
			 
			for(int i = 0; i < other.na; ++i)
				new_atoms[i] = other.a[i];

            // 2: освобождаем "старую" память
            delete [] a;
 
            // 3: присваиваем значения в "новой" памяти объекту
            a = new_atoms;
            na = other.na;

			num_sort = other.num_sort;
			state    = other.state;
			idp      = other.idp;
			nvertex  = other.nvertex;

			L        = other.L;

			topoly = other.topoly;
			fromvertex = other.fromvertex;
			tovertex = other.tovertex;


			rfactor  = other.rfactor;
			point    = other.point;

			xe = other.xe;
			ye = other.ye;
			ze = other.ze;

			r0 = other.r0;

			qtrn = other.qtrn;
			qtrn_r0 = other.qtrn_r0;
			qtrn_norm = other.qtrn_norm;
        }
        // по соглашению всегда возвращаем *this
        return *this;
    }

	Polyhedron(int probe)
	{
		Polyhedron(0.0f, LENEDGE, 0);
	};

	Polyhedron(Vector ri=0.0f, float len = LENEDGE, int probe = 1)  //по умолчанию создаем тетраэдр Si в центре O в вершинах Si - 0, O - 1;
	{
		point.tag = this;
		state2 = -1;

		xe = Vector(1, 0, 0); // ед вектора осей x,y,z внутри тетраэдра
		ye = Vector(0, 1, 0);
		ze = Vector(0, 0, 1);

		if(probe == 1) //1-реальный тетраэдр, 0-пробный
		{
			idp = npall;
			++npall;
		}
		
		r0 = ri;
		state = -1;                      
		nvertex = 4;
		na = 5;
		L = len;						//длина ребра		
		
		float l2 = L * 0.5f;            //половина ребра L=r(O-O)
		float h =  L * sqrtf(2.0f/3);   //высота тетраэдра
		float h1 = L / sqrtf(3.0f);
		float h2 = L * sqrtf(3.0f) / 6;        
		float rt = L * sqrtf(6.0f) / 4; //радиус описаной окружности = r(Si-O)
		float h3 = h - rt;

		a = new Atom[na];

		a[0] = Atom( 0,   0,   0,  Si, 0);
		a[1] = Atom( h1,  0,  -h3, O, 1);
		a[2] = Atom(-h2, -l2, -h3, O, 2);
		a[3] = Atom(-h2,  l2, -h3, O, 3);
		a[4] = Atom( 0,   0,   rt, O, 4);


		for(int i = 0; i < na; ++i)
		{
			a[i].allowtoconnect = 1;
			a[i].connected = 0;
			a[i].deleted = 0;
		}

		a[0].allowtoconnect = 0;
		a[0].addbond(&a[1]);
		a[0].addbond(&a[2]);
		a[0].addbond(&a[3]);
		a[0].addbond(&a[4]);

		a[1].addbond(&a[0]);
		a[1].addbond(&a[2]);
		a[1].addbond(&a[3]);
		a[1].addbond(&a[4]);
		
		a[2].addbond(&a[0]);
		a[2].addbond(&a[1]);
		a[2].addbond(&a[3]);
		a[2].addbond(&a[4]);

		a[3].addbond(&a[0]);
		a[3].addbond(&a[1]);
		a[3].addbond(&a[2]);
		a[3].addbond(&a[4]);

		a[4].addbond(&a[0]);
		a[4].addbond(&a[1]);
		a[4].addbond(&a[2]);
		a[4].addbond(&a[3]);
	}

	//получаем координаты атома i в глоб. системе координат
	Vector lc(int i){ if(i >= 0 && i < na) return a[i].lc; else return -1; }
	
	//получаем координаты атома i в лок. системе координат
	Vector gc(int i){ if(i >= 0 && i < na) return a[i].gc; else return -1; }

	//радиус вектор, установка равносильна перемещению многоранника в точку r0
	inline void set_r0(Vector r)
	{
		r0 = r - r0;
		xe = xe + r0;
		ye = ye + r0;
		ze = ze + r0;

		for(int i = 0; i < na; ++i)
			a[i].gc = a[i].gc + r0;
	}

	//перемещаем многоранник из пар. оси from-to
	inline void movefromto(Vector from, Vector to)
	{
		set_r0(r0 + to - from);
	}

	//присоединяет вершину va этого  тетраэдра к вершине другого коаксиально и паралельно
	inline int connectcoaxandparallel(int va, Polyhedron &toPoly, int vb)
	{
		
		if( a[va].connected == 1 || toPoly.a[vb].connected == 1 || 
			a[va].allowtoconnect == 0 || toPoly.a[vb].allowtoconnect == 0 )
			return -1;

		//random_rotvertex(-15*PI/180, 15*PI/180);  //случ. образом вращаем все вершины тетраэдра вокруг центра

		topoly = &toPoly;
		fromvertex = va;
		tovertex = vb;

		Vector vertexa(a[va].gc);
		Vector vertexb(toPoly.a[vb].gc);
		
		movefromto(vertexa, vertexb);        //совмещаем вершины
		
		vertexa = a[va].gc - r0;
		vertexb = vertexb - toPoly.r0;

		float angle;                         //угол м.у. векторами

		Vector crossab = vertexb * vertexa;  //векторное произведение

		if(crossab.norm3() < TINY)			 //т.е. оси совпадают
		{
			//проверяем  угол, 0-векторы совпадают 180-противопол направлены
			angle = ((vertexa & vertexb ) < 0.0f) ? 0 : PI; 
			
			//создаем препендикулярную ось вручную
			crossab = perpaxisfree(vertexa);

			//проверка на нулевйю ось
			if(crossab.x < TINY && crossab.y < TINY && crossab.x < TINY)
				return -1;
		}
		else
		{
			 angle = PI - anglebetweentwovectors(vertexa, vertexb);
		}

		float adeg = angle * 180 / PI;

		rotate(vertexa + r0, crossab,  angle);

		if(va == 1 || va == 2 || va == 3) //выбираем оси для совмещения
			vertexa = ze - r0;           
		else 
			vertexa = xe - r0;

		if(vb == 1 || vb == 2 || vb == 3) 
			vertexb = toPoly.ze - toPoly.r0;
		else
			vertexb = toPoly.xe - toPoly.r0;

		angle = angvectprojections(a[va].gc - r0, vertexa, vertexb);
		
		rotate(r0, a[va].gc - r0 , angle);	//разворачиваем пристраиваемый тетраэдр чтобы
											//чтобы оси X-X или X-Z совпадали

		a[va].connected = 1;
		a[va].deleted = 1;
		//a[va].replacebonds(&toPoly.a[vb]);

		toPoly.a[vb].connected = 1;

		return 0;
	}

	//вычисление угла м.у. проекциями векторов a и b на плоскость plane
	inline float angvectprojections(const Vector & plane, const Vector & a, const Vector & b)
	{
		float planenorm_inv = 1.0f / plane.norm3();

		float ta = - (plane & a) * planenorm_inv;
		float tb = - (plane & b) * planenorm_inv;

		Vector proja = (plane * ta + a);
		Vector projb = (plane * tb + b);

		float angle = anglebetweentwovectors(proja, projb);

		angle = ( (plane & (proja * projb)) < 0.0f) ? -angle : angle; //знак поворот в какую сторону

		return angle;
	}

	//получаем вектор перепендикулярный к вектору v
	inline Vector perpaxisfree(const Vector & v)
	{
		if (v.x > TINY) 
			return Vector((-v.y - v.z) / v.x, 1.0f, 1.0f);
		else if (v.y > TINY) 
			return  Vector(1.0f, (-v.x - v.z) / v.y, 1.0f);
		else if (v.z > TINY)
			return Vector(1.0f, 1.0f, (-v.x - v.y) / v.z);
		else 
			return Vector(0.0f, 0.0f, 0.0f);
	}
	

	inline int connecttofreevertex(Polyhedron topolyhedron, int & tovertex, int & fromvertex);

	//устраняем коллизии атомов (для глобального массива tetra)
	inline void delcol();

	//устраняем коллизии атомов
	inline void delcol(const Polyhedron *tets, int ntets);
	/*
	void setqtrn(const Vector & r0,  const Vector & axes, const float angle)
	{
		float half_angle = angle*0.5f;
        float sin_half_a = sinf(half_angle); 
		Vector v(axes);
		v.normalize3();
		qtrn_r0 = Vector(r0);
		qtrn = Vector(v.x*sin_half_a, v.y*sin_half_a, v.z*sin_half_a, cosf(half_angle));
		qtrn_norm = qtrn.norm4();
	}
	*/

	//ващение тетраэдра, r - начало оси вращения, ось, угол
	void rotate(Vector r, Vector axes, float angle)
	{
		float half_angle = angle*0.5f;
        float sin_half_a = sinf(half_angle); 

		Vector v(axes);
		v.normalize3();

		qtrn_r0 = r;
		qtrn = Vector(v.x*sin_half_a, v.y*sin_half_a, v.z*sin_half_a, cosf(half_angle));
		
		float qtrn_norm = qtrn.norm4();

		for(int i = 0; i < na; ++i)
		{
			a[i].gc = a[i].gc - qtrn_r0;
			a[i].gc.rotate(qtrn, qtrn_norm);
			a[i].gc = a[i].gc + qtrn_r0;
		}

		xe = xe - qtrn_r0;
		ye = ye - qtrn_r0;
		ze = ze - qtrn_r0;
		xe.rotate(qtrn, qtrn_norm);
		ye.rotate(qtrn, qtrn_norm);
		ze.rotate(qtrn, qtrn_norm);
		xe = xe + qtrn_r0;
		ye = ye + qtrn_r0;
		ze = ze + qtrn_r0;

		r0 = r0 - qtrn_r0;
		r0.rotate(qtrn, qtrn_norm);
		r0 = r0 + qtrn_r0;
		

	}


	inline void rotate_lcs(int num_a, Vector r0_lcs, Vector axes, float angle)
	{
		float half_angle = angle*0.5f;
        float sin_half_a = sinf(half_angle); 

		Vector v(axes);
		v.normalize3();

		Vector qtrn_in(v.x*sin_half_a, v.y*sin_half_a, v.z*sin_half_a, cosf(half_angle));
		
		float norm = qtrn_in.norm4();

		a[num_a].lc = a[num_a].lc - r0_lcs;
		a[num_a].lc.rotate(qtrn_in, norm);
		a[num_a].lc = a[num_a].lc + r0_lcs;
		
		a[num_a].gc = a[num_a].lc + r0;
	}
	
	inline void movetopos_lcs(int num_a, Vector new_pos)
	{
		a[num_a].lc = new_pos;
		a[num_a].gc = a[num_a].lc + r0;
	}

	inline void movebydelta_lcs(int num_a, Vector deltar)
	{
		a[num_a].lc = a[num_a].lc + deltar;
		a[num_a].gc = a[num_a].lc + r0;
	}

	inline void rotatevertex(int vertex_n, Vector axis_r0, Vector axis, float angle)
	{
		float half_angle = angle*0.5f;
        float sin_half_a = sinf(half_angle); 

		Vector v(axis);
		v.normalize3();

		Vector qtrn_in(v.x*sin_half_a, v.y*sin_half_a, v.z*sin_half_a, cosf(half_angle));
		
		float norm = qtrn_in.norm4();

		a[vertex_n].gc = a[vertex_n].gc - axis_r0;
		a[vertex_n].gc.rotate(qtrn_in, norm);
		a[vertex_n].gc = a[vertex_n].gc + axis_r0;
		
		a[vertex_n].lc = a[vertex_n].gc - r0;
	}

	//задаем угол между двумя атомами, ось вращения проходит через центр(r0) многогранника
	inline void setanglebetweentwovertex(int vertexa_n, int vertexb_n, float angle)
	{
		Vector perpaxis = (a[vertexa_n].gc - r0) * (a[vertexb_n].gc - r0);
		
		float angle_delta = (angle - anglebetweentwovectors(a[vertexa_n].gc - r0, a[vertexb_n].gc - r0)) / 2;
		
		//rotate_lcs(vertexa_n, 0, perpaxis, -angle_delta);
		//rotate_lcs(vertexb_n, 0, perpaxis, angle_delta);

		rotatevertex(vertexa_n, 0, perpaxis, -angle_delta);
		rotatevertex(vertexb_n, 0, perpaxis, angle_delta);



	}
		

	inline void randomrotatevertexbysolidangle(int vertex_n, Vector axis_r0, Vector axis, float anglefrom, float angleto)
	{
		float randomangle;
		
		Vector perp_axis = perpaxisfree(axis);

		//вращение перп. оси на произв. угол
		randomangle = rand32fi(0, PI * 2);  
		float half_angle = randomangle * 0.5f;
        float sin_half_a = sinf(half_angle); 
		Vector v(axis);
		v.normalize3();
		qtrn = Vector(v.x*sin_half_a, v.y*sin_half_a, v.z*sin_half_a, cosf(half_angle));
		perp_axis.rotate(qtrn, qtrn.norm4());

		//вращение на случ. угол в пределах телесного угла
		randomangle = rand32fi(anglefrom, angleto); 
		rotatevertex(vertex_n, axis_r0, perp_axis, randomangle);
	}

	inline void randomrotatevertexbysolidangle2(int vertex_n, float anglefrom, float angleto)
	{
		float randomangle;
		Vector axis = a[vertex_n].gc - r0;
		Vector perp_axis = perpaxisfree(axis);

		//вращение перп. оси на произв. угол
		randomangle = rand32fi(0, PI * 2);  
		float half_angle = randomangle * 0.5f;
        float sin_half_a = sinf(half_angle); 
		Vector v(axis);
		v.normalize3();
		qtrn = Vector(v.x*sin_half_a, v.y*sin_half_a, v.z*sin_half_a, cosf(half_angle));
		perp_axis.rotate(qtrn, qtrn.norm4());

		//вращение на случ. угол в пределах телесного угла
		randomangle = rand32fi(anglefrom, angleto); 
		rotatevertex(vertex_n, r0, perp_axis, randomangle);
	}



}; //class Polehedron


extern int Ntetra;
extern Polyhedron *tetra;


inline void Polyhedron :: delcol()
{
	float r, minR;
	int sa1, sa2;
	
	Polyhedron &testp = *this;

	for(int ip = 0; ip < testp.na; ++ip)
	{
		if(testp.a[ip].deleted == 1) continue;

		for(int it = 0; it < Ntetra; ++it)
		{
			for(int ita = 0; ita < tetra[it].na; ++ita)
			{
				if(tetra[it].a[ita].deleted == 1)continue;

				r = dr(tetra[it].a[ita].gc, testp.a[ip].gc);
				
				sa1 = testp.a[ip].sort;
				sa2 = tetra[it].a[ita].sort;

				minR = atomscollen[sa1][sa2];

				if(r < minR)
				{
					testp.a[ip].deleted = 1;
					//testp.a[ip].replacebonds(&(tetra[it].a[ita]));
				}

			}
		}
	}
}



inline void Polyhedron :: delcol(const Polyhedron *tets, int ntets)
{
	float r, minR;
	int sa1, sa2;
	
	Polyhedron &testp = *this;

	for(int i = 0; i < testp.na; ++i)
	{
		if(testp.a[i].deleted == 1) continue;

		for(int j = 0; j < ntets; ++j)
		{
			if(tets[j].state2 == -1) continue;

			for(int k = 0; k < tets[j].na; ++k)
			{
				if(tets[j].a[k].deleted == 1)continue;

				r = dr(tets[j].a[k].gc, testp.a[i].gc);
				
				sa1 =   testp.a[i].sort;
				sa2 = tets[j].a[k].sort;

				minR = atomscollen[sa1][sa2];

				if(r < minR)
				{
					testp.a[i].deleted = 1;
					//testp.a[i].replacebonds(&(tets[j].a[k]));
				}

			}
		}
	}
}

inline int Polyhedron :: connecttofreevertex(Polyhedron topolyhedron, int & tovertex, int & fromvertex)
	{
		for(int i = 0; i < na; ++i)
		{
			if(a[i].connected == 0 && a[i].allowtoconnect == 1)
			{
				for(int j = 0; j < topolyhedron.na; ++j)
				{
					if(topolyhedron.a[j].connected == 0 && topolyhedron.a[j].allowtoconnect == 1)
					{
						if(0 == connectcoaxandparallel(i, topolyhedron, j))
						{
							delcol();
							fromvertex = i;
							tovertex = j;
							return 0;
						}
					}
				}
			}
		}
	
		return -1;
	}



inline void writetetra(string filename)
{
	ofstream outf(filename);
	
	int row = 1;

	for(int p = 0; p < Ntetra; ++p)
		for(int i = 0; i < tetra[0].na; ++i)
		{
			if(tetra[p].a[i].deleted == 1) continue;

			
			tetra[p].a[i].ida = row;
			outf << setw(5) << tetra[p].a[i].anum;
			outf << setw(25) << setprecision(4) << fixed;
			outf << tetra[p].a[i].gc.x;
			outf << setw(25) << setprecision(4) << fixed;
		    outf << tetra[p].a[i].gc.y;
			outf << setw(25) << setprecision(4) << fixed;
			outf << tetra[p].a[i].gc.z << endl;

			++row;
		}
	

}

inline void writetetra2(string filename, const Polyhedron &poly)
{
	ofstream outf(filename);
	
	int row = 1;

	tetra[Ntetra] = poly;
	++Ntetra;


	for(int p = 0; p < Ntetra; ++p)
		for(int i = 0; i < tetra[0].na; ++i)
		{
			if(tetra[p].a[i].deleted == 1) continue;

			
			tetra[p].a[i].ida = row;
			outf << setw(5) << tetra[p].a[i].anum;
			outf << setw(25) << setprecision(4) << fixed;
			outf << tetra[p].a[i].gc.x;
			outf << setw(25) << setprecision(4) << fixed;
		    outf << tetra[p].a[i].gc.y;
			outf << setw(25) << setprecision(4) << fixed;
			outf << tetra[p].a[i].gc.z << endl;

			++row;
		}
	
		--Ntetra;

}

inline void writetetra3(string filename, const Polyhedron *tets, int ntets)
{
	ofstream outf(filename);
	
	int row = 1;

	for(int p = 0; p < ntets; ++p)
		for(int i = 0; i < tets[0].na; ++i)
		{
			if(tets[p].a[i].deleted == 1) continue;
			
			tets[p].a[i].ida = row;
			outf << setw(5) << tets[p].a[i].anum;
			outf << setw(25) << setprecision(4) << fixed;
			outf << tets[p].a[i].gc.x;
			outf << setw(25) << setprecision(4) << fixed;
		    outf << tets[p].a[i].gc.y;
			outf << setw(25) << setprecision(4) << fixed;
			outf << tets[p].a[i].gc.z << endl;

			++row;
		}
	
}