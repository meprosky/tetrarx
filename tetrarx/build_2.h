#pragma once

#include"nvector.h"
#include"param.h"
#include"polyhedron.h"
#include"optclass.h"



class OnePoly
{
public:
 
	int ntets;
	Polyhedron *tets;

	DXRCL dxrcl_cluster;

	OPTIMIZE rosenbrock;
	int t0, v0, v1;
	int pvtdim;
	nVector pvt, pvtmin, pvtmax, dpvt, dpvtstep, pvtminrosenbrock, pvtmaxrosenbrock; 
	Polyhedron poly;

	OnePoly() 
	{
	}

	void set_dxrcl(string filename)
	{
		dxrcl_cluster.dxrclinit(filename);
	}

	void set_tets(Polyhedron *tets_, int ntets_)
	{
		if(ntets > 0 && tets != nullptr)  delete [] tets;
		ntets = ntets_ + 1;
		tets = new Polyhedron[ntets];


		for(int i = 0; i < ntets-1; ++i)
		{
			tets[i].state2 = 1;
			tets[i] = tets_[i];
		}

	}

	void set_pvt(int pvtdim_, nVector pvtmin_, nVector pvtmax_, float step1, float step2)
	{ 
		pvtdim = pvtdim_; pvtmin = pvtmin_; pvtmax = pvtmax_; 
		pvt.initnew(2);
		dpvt = pvtmax - pvtmin;

		dpvtstep.initnew(pvtdim);

		dpvtstep.x[0] = step1;
		dpvtstep.x[1] = step2;

		pvtminrosenbrock = -pvtmax;  //из-за цикличности углов
		pvtmaxrosenbrock = pvtmax;
	}

	inline void set_rosenbrock_minmax(const nVector &min_, const nVector &max_)
	{
		pvtminrosenbrock = min_;
		pvtmaxrosenbrock = max_;
	}

	void set_rosenbrockparameters(
		float epsfunc_, float epssteps_, 
		float step_, float k_alpha_, float k_beta_,
		int isusepenalty_,
		const int maxitertation_,
		float (*func_)(const nVector &p, void *tag))
	{
		nVector steps(8);
		steps = step_;  // пока шаг один для всех

		rosenbrock.setparameters(pvtdim, epsfunc_, epssteps_, steps, k_alpha_, k_beta_, 
			pvtminrosenbrock, pvtmaxrosenbrock, isusepenalty_, maxitertation_, func_);
	}


	//выполнить процедуру поиска минимума по Розенброку
	void execrosenbrock(nVector &pointopt, float &fvalopt)
	{
		pvt.tag = this;
		rosenbrock.findoptimum(pvt, pointopt, fvalopt);
	}

	

	inline void set_randpvtvector_andstep()
	{

		for(int i = 0; i < pvtdim; ++i)
		{
			int ns = (int)(dpvt.x[i] / dpvtstep.x[i]) + 1;
			int rnds  = rand32ui(0, ns);

			pvt.x[i] = rand32fi(pvtmin.x[i], pvtmin.x[i] + dpvtstep.x[i] * rnds);
		}
	}

	void rand_attempt(int nattempts, float &fvalopt_)
	{
		nVector point(2), pointopt(2), pointcur(2);
		float fvalopt = 1000.0f, fvalcur;
		int topt=0, vopt=0;

		for(int i = 0; i < nattempts; ++i)
		{
			Polyhedron probe;
			int t = rand32ui(0, ntets - 2);
			int v = rand32ui(1, 4); 
			//int v1 = rand32ui(1, 4); 
			
			set_randpvtvector_andstep();
			
			probe.state2 = 1;
			probe.connectcoaxandparallel(1, tets[t], v);
			probe.rotate(probe.gc(1), probe.gc(0) - probe.gc(1), pvt.x[0]);
			probe.rotate(probe.gc(1), probe.ye    - probe.r0,    pvt.x[1]);
			probe.delcol(tets, ntets-2);

			tets[ntets-1] = probe;
			tets[ntets-1].state2 = 1;

			writetetra3("__temp.xyz", tets, ntets);

			fvalcur = dxrcl_cluster.rfactor_wcalcdxrcl(tets, ntets);

			if(fvalcur < fvalopt)
			{
				fvalopt = fvalcur;
				pointopt = pointcur;
				t0 = t;
				v0 = v;
				v1 = 1;
			}

		}

		writetetra3("__temp.xyz", tets, ntets);

		execrosenbrock(pointcur, fvalcur);

		writetetra3("__temp.xyz", tets, ntets);

	}


};


float func_forrosenbrock2(const nVector &point, void *tag)
{
	OnePoly &op = *(static_cast<OnePoly*>(tag));
	
	Polyhedron probe;

	probe.state2 = 1;
	probe.connectcoaxandparallel(op.v1, op.tets[op.t0], op.v0);
	probe.rotate(probe.gc(1), probe.gc(0) - probe.gc(1), point.x[0]);
	probe.rotate(probe.gc(1), probe.ye    - probe.r0,    point.x[1]);
	probe.delcol(op.tets, op.ntets-2);

	op.tets[op.ntets-1] = probe;

	float rff = op.dxrcl_cluster.rfactor_wcalcdxrcl(op.tets, op.ntets);

	return rff; 
}




class PolyhedronMerge
{
public:
	
	
	
	int n0;
	Polyhedron *tets0;
	
	int n1;
	Polyhedron *tets1;

	int nm;
	Polyhedron *tets;

	int t0, v0, t1, v1;

	PolyhedronMerge()
	{
		n0 = 0; tets0 = nullptr;
		n1 = 0; tets1 = nullptr;
		nm = 0; tets = nullptr;
		t0 = -1; v0 = -1; t1 = -1; v1 = -1;
	}

	void reset_tets0(int n0_)
	{
		if(n0 > 0 && tets0 != nullptr)  delete [] tets0;
		n0 = n0_;
		tets0 = new Polyhedron[n0];
	}

	void reset_tets1(int n1_)
	{
		if(n1 > 0 && tets1 != nullptr)  delete [] tets1;
		n1 = n1_;
		tets1 = new Polyhedron[n1];
	}

	void reset_tets(int nm_)
	{
		if(nm > 0 && tets != nullptr)  delete [] tets;
		nm = nm_;
		tets = new Polyhedron[nm];
	}

	void set_t0(Polyhedron *tets0_, int n0_)
	{
		reset_tets0(n0_);
		for(int i = 0; i < n0; ++i) 	tets0[i] = tets0_[i];
	}

	void set_t1(Polyhedron *tets1_, int n1_)
	{
		reset_tets1(n1_);
		for(int i = 0; i < n1; ++i) 	tets1[i] = tets1_[i];
	}

	void merge_tets()
	{
		reset_tets(n0+n1);
		for(int i = 0;  i < n0; ++i) 
			tets[i] = tets0[i];
		
		for(int i = n0; i < nm; ++i) 
			tets[i] = tets1[i-n0];
	}

	void merge_tets_to_t0()
	{
		merge_tets();
		if(n0 > 0 && tets0 != nullptr)  delete [] tets0;
		tets0 = new Polyhedron[nm];
		for(int i = 0; i < nm; ++i) tets0[i] = tets[i];
	}

	void connect_tets1to0(int from_t1, int from_v1, int to_t0, int to_v0)
	{
		t0 = to_t0; v0 = to_v0; t1 = from_t1; v1 = from_v1;

		Vector p0 = tets0[t0].a[v0].gc;
		Vector p1 = tets1[t1].a[v1].gc;

		reset_tets(n0+n1);

		for(int i = 0;  i < n0; ++i) 
		{
			tets[i].state2 = 1;
			tets[i] = tets0[i];
		}
		
		int j;
		for(int i = n0;  i < nm; ++i)
		{
			j = i - n0;
			tets1[j].state2 = 1;
			tets1[j].movefromto(p1, p0);
			tets1[j].delcol(tets, i);
			 
			tets[i] = tets1[j];
		}
	}

	void move_tets1(Vector from, Vector to)
	{
		reset_tets1(n1);
		for(int i = 0;  i < n1; ++i) { tets1[i].state2 = 1; tets1[i].movefromto(from, to);}
	}

	void rotate_tets0(const Vector &r0, const Vector &axis, float angle)
	{
		reset_tets0(n0);
		for(int i = 0;  i < n0; ++i){ tets0[i].state2 = 1; tets0[i].rotate(r0, axis, angle);}
	}

	void rotate_tets1(const Vector &r0, const Vector &axis, float angle)
	{
		reset_tets1(n1);
		for(int i = 0;  i < n1; ++i){ tets1[i].state2 = 1; tets1[i].rotate(r0, axis, angle);}
	}

};


