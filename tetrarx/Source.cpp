#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<time.h>
#include<string>
#include<vector>
#include"polyhedron.h"
#include"param.h"
#include"dxrcl.h"
#include"optclass.h"
#include"rfactoropt.h"
#include"mergesort.h"
#include"build_1.h"
#include"build_2.h"

using namespace std;

//глобальные перменные
int Ntetra;
Polyhedron *tetra = new Polyhedron[NUMSTETRA];
Polyhedron *temp  = new Polyhedron[PROBETETRA];

DXRCL dxrcl; 

int Polyhedron::npall = 0;

//void _formclust1();
//void _formclust2();
void _formclust3();

void main()
{
	setlocale(LC_CTYPE, "");			//чтобы в консоли были русские буквы

	DXRCL ldxrcl("__4t_sio2_hs.dxr"); ;
	dxrcl = ldxrcl;

	atomscollen[Si][O]  = 1.4f;
	atomscollen[O][Si]  = 1.4f;
	atomscollen[O][O]   = 2.4f;
	atomscollen[Si][Si] = 2.2f;

	/*
	time_t t;
	t = clock();
	//...
	t = clock() - t;
	cout << "t1 = " << t  << endl;
	t = clock();
	//...
	t = clock() - t;
	cout << "t2 = " << t << endl;
	*/
	
	
	tetra[0].state = 1;            //состояние - активный 
	tetra[0].set_r0(Vector(0,0,0)); //помещаем в начало координат, хотя по умолчанию все в начале координат
	Ntetra = 1;                    // число тетраэдров в системе
	

	//formclust2(); //нова конфигурация старые 4 + 2 новых тетраэдра
	_formclust3();
	exit(1);

	nVector pmin(2), pmax(2);
	pmin.x[0] = 0.0f  *PI/180.0f;   pmin.x[1] = 0.0f*PI/180.0f;
	pmax.x[0] = 360.0f*PI/180.0f;   pmax.x[1] = 45.0f*PI/180.0f;
		
	BuildClustConf bclust;
	//bclust.set_dxrcl("__4t_sio2_hs.dxr");
	bclust.set_dxrcl("__6t_sio2_hs.dxr");
	bclust._exp1();
	bclust.set_pvtdimpertetra(2);
	bclust.set_pvtminmax_oneforall(pmin, pmax);
	

	float fvalcalc, fvalopt = 10000.0f;
	nVector popt(8);

	bclust.set_dstep_2dim(5.0f*PI/180.0f, 5.0f*PI/180.0f); 

	for(int i = 0; i < 200000; ++i)
	{
		bclust.build_wrandcluster_wpvtminmaxoneforall_wsteps();
		fvalcalc = bclust.rfactorcalc();

		if(fvalcalc < fvalopt)
		{
			fvalopt = fvalcalc;
			popt = bclust.pvt;
		}
	}

    bclust.build_cluster_wpvtvector(popt);
	fvalopt = bclust.rfactorcalc();

	//bclust.build_wrandcluster_wpvtminmaxoneforall();
	
	nVector prosen(8);
	float fvalrosen;


	bclust.set_rosenbrockparameters(0.1e-8f, 0.1e-8f, 0.1f, 3.0f, -0.5f, 1, 500, func_forrosenbrock);
	bclust.rosenbrock.setiterationmode(1);
	bclust.execrosenbrock(prosen, fvalrosen);

	bclust.write_tetstofile("__temp.xyz");
	bclust.dxrcl_cluster.Write_Hs("__temphscalc.xy", bclust.dxrcl_cluster.hs_calc_normby1);
	bclust.dxrcl_cluster.Write_Hs("__temphsexp.xy", bclust.dxrcl_cluster.hs_exp_normby1);


	for(int i = 0; i < bclust.ntets; ++i)
	{
		Polyhedron &poly = bclust.tets[i];
		for(int j = 0; j < poly.na; ++j)
		{
			if(poly.a[j].allowtoconnect == 0 || poly.a[j].deleted == 1 || poly.a[j].connected == 1)
				continue;


			//bclust.add(i, j, bclust.ntets

		}


	}




	exit(1);

	

	float optrf = 1000.0f;
	nVector optpoint(8);



	for(int i = 0; i < 15000; ++i)
	{
		bclust.build_wrandcluster_wpvtminmaxoneforall();
		
		nVector prosen(8);
		float fvalrosen;
		//bclust.execrosenbrock(prosen, fvalrosen);

		fvalrosen = bclust.rfactorcalc();

		if(fvalrosen < optrf)
		{
			optrf = fvalrosen;
			//optpoint = prosen;
			optpoint = bclust.pvt;
		}
	}

	//nVector prosen(8);
	//float fvalrosen;
	//bclust.execrosenbrock(prosen, fvalrosen);


	//bclust.build_cluster_wpvtvector(optpoint);
	bclust.write_tetstofile("__temp.xyz");
	bclust.dxrcl_cluster.Write_Hs("__temphscalc.xy", bclust.dxrcl_cluster.hs_calc_normby1);
	bclust.dxrcl_cluster.Write_Hs("__temphsexp.xy", bclust.dxrcl_cluster.hs_exp_normby1);

	exit(1);

	delete [] tetra;

}









void _formclust1()
{
	Ntetra = 6;
	tetra[1].connectcoaxandparallel(1, tetra[0], 1);
	tetra[1].rotate(tetra[1].gc(1), tetra[1].gc(0) - tetra[1].gc(1), PI*60/180);
	tetra[1].rotate(tetra[1].gc(1), tetra[1].ye - tetra[1].r0, PI*30/180);

	tetra[2].connectcoaxandparallel(1, tetra[0], 2);
	tetra[2].rotate(tetra[2].gc(1), tetra[2].gc(0) - tetra[2].gc(1), PI*60/180);
	tetra[2].rotate(tetra[2].gc(1), tetra[2].ye - tetra[2].r0, PI*30/180);

	tetra[3].connectcoaxandparallel(1, tetra[1], 4);
	tetra[3].rotate(tetra[3].gc(1), tetra[3].gc(0) - tetra[3].gc(1), PI*60/180);
	tetra[3].rotate(tetra[3].gc(1), tetra[3].ye - tetra[3].r0, PI*30/180);
	
	//еще плюс 2
	tetra[4].connectcoaxandparallel(1, tetra[3], 4);
	tetra[4].rotate(tetra[4].gc(1), tetra[4].gc(0) - tetra[4].gc(1), PI*90/180);
	tetra[4].rotate(tetra[4].gc(1), tetra[4].ye    - tetra[4].r0, PI*45/180);

	tetra[5].connectcoaxandparallel(1, tetra[4], 4);
	tetra[5].rotate(tetra[5].gc(1), tetra[5].gc(0) - tetra[5].gc(1), PI*90/180);
	tetra[5].rotate(tetra[5].gc(1), tetra[5].ye    - tetra[5].r0, PI*45/180);


	writetetra("__newtemp.xyz");
	dxrcl.dxrclcalc();



	dxrcl.Write_Hs("__newtemp6hscalc.xy", dxrcl.hs_calc);
	dxrcl.Write_Hs("__newtemp6hsnorm.xy", dxrcl.hs_calc_normby1);


};

/*
void _formclust2()
{
	int n1 = 2, n2=3;
	Polyhedron *t1 = new Polyhedron[2];
	Polyhedron *t2 = new Polyhedron[3];

	t1[1].connectcoaxandparallel(1, t1[0], 1);
	t1[1].rotate(t1[1].gc(1), t1[1].gc(0) - t1[1].gc(1), PI*60/180);
	t1[1].rotate(t1[1].gc(1), t1[1].ye - t1[1].r0, PI*30/180);

	t2[0].set_r0(Vector(10, 10, 10));

	t2[1].connectcoaxandparallel(1, t2[0], 1);
	t2[1].rotate(t2[1].gc(1), t2[1].gc(0) - t2[1].gc(1), PI*60/180);
	t2[1].rotate(t2[1].gc(1), t2[1].ye - t2[1].r0,       PI*30/180);

	t2[2].connectcoaxandparallel(1, t2[1], 4);
	t2[2].rotate(t2[2].gc(1), t2[2].gc(0) - t2[2].gc(1), PI*90/180);
	t2[2].rotate(t2[2].gc(1), t2[2].ye    - t2[2].r0,    PI*45/180);

	PolyhedronMerge p;

	p.set_t0(t1, n1);
	p.set_t1(t2, n2);
	p.merge_tets();
	p.connect_tets1to0(0, 4, 1, 3);
	//p.rotate_tets1(p.tets1[0].gc(4), p.tets1[0].ye - p.tets1[0].r0, -PI*90/180); 

	writetetra3("__temp2.xyz", p.tets, p.nm);
};
*/

void _formclust3()
{
	Polyhedron *t1 = new Polyhedron[4];
	t1[1].connectcoaxandparallel(1, t1[0], 1);
	t1[1].rotate(t1[1].gc(1), t1[1].gc(0) - t1[1].gc(1), PI*60/180);
	t1[1].rotate(t1[1].gc(1), t1[1].ye - t1[1].r0, PI*30/180);

	t1[2].connectcoaxandparallel(1, t1[0], 2);
	t1[2].rotate(t1[2].gc(1), t1[2].gc(0) - t1[2].gc(1), PI*60/180);
	t1[2].rotate(t1[2].gc(1), t1[2].ye - t1[2].r0, PI*30/180);

	t1[3].connectcoaxandparallel(1, t1[1], 4);
	t1[3].rotate(t1[3].gc(1), t1[3].gc(0) - t1[3].gc(1), PI*60/180);
	t1[3].rotate(t1[3].gc(1), t1[3].ye - t1[3].r0, PI*30/180);

	OnePoly op;

	nVector pmin(2), pmax(2);
	pmin.x[0] = 0.0f  *PI/180.0f;   pmin.x[1] = 0.0f*PI/180.0f;
	pmax.x[0] = 360.0f*PI/180.0f;   pmax.x[1] = 45.0f*PI/180.0f;

	op.set_pvt(2, pmin, pmax, 5.0f*PI/180.0f, 5.0f*PI/180.0f);
	op.set_tets(t1, 4);
	op.set_dxrcl("__6t_sio2_hs.dxr");
	op.set_rosenbrockparameters(0.1e-8f, 0.1e-8f, 0.1f, 3.0f, -0.5f, 1, 500, func_forrosenbrock2);
	op.rosenbrock.setiterationmode(1);

	float fvalopt;
	op.rand_attempt(200, fvalopt);



}


// тестовые 4 тетраэдра по которым проверяем алгоритм
    /*
	//Ntetra = 4;
	tetra[1].connectcoaxandparallel(1, tetra[0], 1);
	tetra[1].rotate(tetra[1].gc(1), tetra[1].gc(0) - tetra[1].gc(1), PI*60/180);
	tetra[1].rotate(tetra[1].gc(1), tetra[1].ye - tetra[1].r0, PI*30/180);

	dxrcl.Write_Hs("__temp2thsnorm.xy", dxrcl.hs_calc_normby1);
	dxrcl.Write_Hs("__temp4thsnorm.xy", dxrcl.hs_exp_normby1);

	tetra[2].connectcoaxandparallel(1, tetra[0], 2);
	tetra[2].rotate(tetra[2].gc(1), tetra[2].gc(0) - tetra[2].gc(1), PI*60/180);
	tetra[2].rotate(tetra[2].gc(1), tetra[2].ye - tetra[2].r0, PI*30/180);

	tetra[3].connectcoaxandparallel(1, tetra[1], 4);
	tetra[3].rotate(tetra[3].gc(1), tetra[3].gc(0) - tetra[3].gc(1), PI*60/180);
	tetra[3].rotate(tetra[3].gc(1), tetra[3].ye - tetra[3].r0, PI*30/180);
	*/
