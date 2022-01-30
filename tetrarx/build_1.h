#pragma once

#include"nvector.h"
#include"param.h"
#include"polyhedron.h"
#include"optclass.h"

class TetraedrVertex
{
public:
	int t0, t1, v0, v1;
};

class BuildClustConf
{
public:
	int nseq;					//число пар присоединений тетраэдра-тетраэдр (1-2, 2-1 ...)
	int ntets;					//число тетраэдров
	TetraedrVertex *tv;			//массив номер тетраедра-номер вершины
	Polyhedron *tets;			//тетраэдры

	float rfactor;				//расчитанный R-фактор кластера
	DXRCL dxrcl_cluster;		//расчет dxrcl

	OPTIMIZE rosenbrock;		//процедура поиска минимума по Розенброку

	int pvtdim;					//количество варьируемых параметров
	int pvtdimpertetra;			//количество варьируемых параметров на тетраэдр
	nVector pvt;				//вектор варьируемых параметров
	nVector pvtmin, pvtmax;		//диапазоны варьируемых параметров, на каждый тетраэдр
	nVector pvtminrosenbrock,   //диапазоны варьируемых параметров, для процедуры Розенброка
		pvtmaxrosenbrock;
	
	nVector dpvt;               //pvtmax - pvtmin
	nVector dpvtstep;           //шаг изменения варьируемых параметров в диапазоне dpvt

	nVector pvtmin_oneforall,	//диапазоны варьируемых параметров одни для всех тетраэров
		pvtmax_oneforall;

	//конструктор
	BuildClustConf() { nseq = 0; pvtdimpertetra = -1; tv = nullptr; tets = nullptr; };

    //временная процедура построения последовательностей класетра
	void _exp1()
	{
		add(0,1,1,1);
		add(0,2,2,1);
		add(1,4,3,1);
	}

	//построение последовательностей класетра из внешнего 2-мерного массива
	void build_sequens(int *asq, int nsq)
	{
		for(int i = 0; i < nsq; ++i)
		{
			int idx = i * nsq;
			add( asq[idx], asq[idx+1], asq[idx+2], asq[idx+3]);
		}
	}

	//добавляем последовательность
	void add(int t0, int v0, int t1, int v1)
	{
		int nnew = nseq + 1;
		
		TetraedrVertex *tvnew = new TetraedrVertex[nnew];

		if(nseq > 0 && tv != nullptr) {
			memcpy(tvnew, tv, nseq * sizeof(TetraedrVertex));
			delete [] tv;
			++ntets;
		}
		else
			ntets = 2;

		tvnew[nnew - 1].t0 = t0;
		tvnew[nnew - 1].v0 = v0;
		tvnew[nnew - 1].t1 = t1;
		tvnew[nnew - 1].v1 = v1;

		tv = tvnew;
		nseq = nnew;
	}

	//получаем номера кластера и вершины t0,v0 - тетраэдр к которому присоединяем
	//t1, v1 - номер тетраэдра и вершиный который присоединяем
	inline void get_conf(int c, int &t0, int &v0, int &t1, int &v1)
	{
		if(c >= 0 && c < nseq)	{ t0 = tv[c].t0; v0 = tv[c].v0; t1 = tv[c].t1; v1 = tv[c].v1; }
		else 				    { t0 = -1; v0 = -1; t1 = -1; v1 = -1; }
	}
	
	//задаем число варьируемых параметров на один тетраэдр
	//расчитываем общще число параметров и выделем память в векторе пераметров
	inline void set_pvtdimpertetra(int n_) 
	{
		pvtdim = ntets * n_; 
		pvtdimpertetra = n_;
		pvt.initnew(pvtdim);
		dpvt.initnew(pvtdim);
	}
	
	//задаем вектор варьируемых параметров
	inline void set_pvt(nVector pvt_) 
	{
		pvtdim = pvt_.n;
		pvt = pvt_;
	}

	//задаем границы диапазонов варьируемых параметров
	inline void set_pvtminmax_oneforall(nVector pvtmin_, nVector pvtmax_)
	{
		pvtmin_oneforall = pvtmin_;
		pvtmax_oneforall = pvtmax_;

		pvtmin.initnew(pvtdim);
		pvtmax.initnew(pvtdim);
		pvtminrosenbrock.initnew(pvtdim);
		pvtmaxrosenbrock.initnew(pvtdim);

		for(int i = 0; i < pvtdim; i += pvtdimpertetra)
		{
			pvtmin.x[i]   = pvtmin_oneforall.x[0];
			pvtmin.x[i+1] = pvtmin_oneforall.x[1];

			pvtmax.x[i]   = pvtmax_oneforall.x[0];
			pvtmax.x[i+1] = pvtmax_oneforall.x[1];

			dpvt.x[i] = pvtmax_oneforall.x[0] - pvtmin_oneforall.x[0];
			dpvt.x[i+1] = pvtmax_oneforall.x[1] - pvtmin_oneforall.x[1];
		}


		pvtminrosenbrock = -pvtmax;
		pvtmaxrosenbrock =  pvtmax;

	}
	
	//шаг изменения врьируемых параметров
	inline void set_dstep_2dim(float s1, float s2)
	{
		dpvtstep.initnew(pvtdim);

		for(int i = 0; i < pvtdim; i += pvtdimpertetra)
		{
			dpvtstep.x[i] = s1;
			dpvtstep.x[i+1] = s2;
		}
	}

	//параметры min-max для процедуры Розенброка
	inline void set_rosenbrock_minmax(const nVector &pmin, const nVector &pmax)
	{
		pvtminrosenbrock = pmin;
		pvtmaxrosenbrock = pmax;
	}

	//случайный вектор параметров в заданном диапазоне(диапазоны индивидуальные для каждого тетраэдров)
	inline void set_randpvtvector()
	{
		for(int i = 0; i < pvtdim; ++i)
		{
			pvt.x[i] = rand32fi(pvtmin.x[i], pvtmax.x[i]);
		}
	}

	//генерация случ. вектора в заданном диапазоне и с шагом
	inline void set_randpvtvector_andstep()
	{

		for(int i = 0; i < pvtdim; ++i)
		{
			int ns = (int)(dpvt.x[i] / dpvtstep.x[i]) + 1;
			int rnds  = rand32ui(0, ns);

			pvt.x[i] = rand32fi(pvtmin.x[i], pvtmin.x[i] + dpvtstep.x[i] * rnds);
		}
	}

	//случайный вектор параметров в заданном диапазоне (диапазоны одни для всех тетраэдров)
	inline void set_randpvtvector_fromminmax_oneforall()
	{
		for(int i = 0; i < pvtdim; i += pvtdimpertetra)
		{
			pvt.x[i]   = rand32fi(pvtmin_oneforall.x[0], pvtmax_oneforall.x[0]);
			pvt.x[i+1] = rand32fi(pvtmin_oneforall.x[1], pvtmax_oneforall.x[1]);
		}
	}

	//строим кластер, оси Si-O соосны, другие Si-O паралельны друг другу
	inline void build_init_cluster()
	{
		if(tets != nullptr) delete [] tets;

		ntets = nseq + 1;
		
		tets = new Polyhedron[ntets];
		nVector point(2);
		int t0, v0, t1, v1;

		pvt = 0.0f;

		tets[0].state2 = 1;

		for(int i = 0; i < nseq; ++i)
		{
			get_conf(i, t0, v0, t1, v1);

			Polyhedron &probe = tets[t1];

			probe.connectcoaxandparallel(v1, tets[t0], v0);
			probe.delcol(tets, ntets);
			
			probe.state2 = 1;
			
			point.x[0] = 0.0f;
			point.x[1] = 0.0f;

			probe.point = point;

		}

	}

	//строим случайный кластер
	inline void build_wrandcluster_wpvtminmaxoneforall()
	{
		if(tets != nullptr) delete [] tets;

		ntets = nseq + 1;

		tets = new Polyhedron[ntets];
		nVector point(2);
		int t0, v0, t1, v1, pvtindex;

		set_randpvtvector();

		tets[0].state2 = 1;

		for(int i = 0; i < nseq; ++i)
		{
			get_conf(i, t0, v0, t1, v1);
			
			Polyhedron &probe = tets[t1];

			probe.connectcoaxandparallel(v1, tets[t0], v0);

			pvtindex = t1 * pvtdimpertetra;

			probe.rotate(probe.gc(1), probe.gc(0) - probe.gc(1), pvt.x[pvtindex]);
			probe.rotate(probe.gc(1), probe.ye    - probe.r0,    pvt.x[pvtindex + 1]);
			probe.delcol(tets, ntets);
			
			probe.state2 = 1;
			
			point.x[0] = pvt.x[pvtindex];
			point.x[1] = pvt.x[pvtindex + 1];

			probe.point = point;
		}
	}
	
	inline void build_wrandcluster_wpvtminmaxoneforall_wsteps()
	{
		set_randpvtvector_andstep();
		build_cluster_wpvtvector();
	}

	//строим кластер по параметрам поля pvt
	inline void build_cluster_wpvtvector()
	{
		if(ntets > 0) delete [] tets;
		
		ntets = nseq + 1;
		
		tets = new Polyhedron[ntets];
		nVector point(2);
		int t0, v0, t1, v1, pvtindex;

		tets[0].state2 = 1;

		for(int i = 0; i < nseq; ++i)
		{
			get_conf(i, t0, v0, t1, v1);

			Polyhedron &probe = tets[t1];

			probe.connectcoaxandparallel(v1, tets[t0], v0);

			pvtindex = t1 * pvtdimpertetra;

			probe.rotate(probe.gc(1), probe.gc(0) - probe.gc(1), pvt.x[pvtindex]);
			probe.rotate(probe.gc(1), probe.ye    - probe.r0,    pvt.x[pvtindex + 1]);
			probe.delcol(tets, ntets);
			probe.state2 = 1;
			
			point.x[0] = pvt.x[pvtindex];
			point.x[1] = pvt.x[pvtindex + 1];

			probe.point = point;
		}
	}

	//строим кластер по параметрам получаемым извне
	inline void build_cluster_wpvtvector(const nVector &pvt_)
	{
		if(ntets > 0) delete [] tets;
		
		ntets = nseq + 1;
		
		tets = new Polyhedron[ntets];
		nVector point(2);
		int t0, v0, t1, v1, pvtindex;

		tets[0].state2 = 1;

		pvt = pvt_;

		for(int i = 0; i < nseq; ++i)
		{
			get_conf(i, t0, v0, t1, v1);

			Polyhedron &probe = tets[t1];

			probe.connectcoaxandparallel(v1, tets[t0], v0);

			pvtindex = t1 * pvtdimpertetra;

			probe.rotate(probe.gc(1), probe.gc(0) - probe.gc(1), pvt_.x[pvtindex]);
			probe.rotate(probe.gc(1), probe.ye    - probe.r0,    pvt_.x[pvtindex + 1]);
			probe.delcol(tets, ntets);
			probe.state2 = 1;
			
			point.x[0] = pvt_.x[pvtindex];
			point.x[1] = pvt_.x[pvtindex + 1];

			probe.point = point;
		}
	}

	//инициализация класса расчета интенсивностей рассения
	void set_dxrcl(string filename)
	{
		dxrcl_cluster.dxrclinit(filename);
	}

	//расчет R-фактора кластера
	float rfactorcalc()
	{
		rfactor = dxrcl_cluster.rfactor_wcalcdxrcl(tets, ntets);
		return rfactor;
	}

	//задаем параметры поиска минимума по Розенброку
	void set_rosenbrockparameters(
		float epsfunc_, float epssteps_, 
		float step_, float k_alpha_, float k_beta_,
		int isusepenalty_,
		const int maxitertation_,
		float (*func_)(const nVector &p, void *tag))
	{
		nVector steps(8);
		steps = step_;

		rosenbrock.setparameters(
			pvtdim, 
			epsfunc_, epssteps_, steps, k_alpha_, k_beta_, 
			pvtminrosenbrock, pvtmaxrosenbrock, isusepenalty_, maxitertation_, func_);

		//cout << "ssss";

	}

	//выполнить процедуру поиска минимума по Розенброку
	void execrosenbrock(nVector &pointopt, float &fvalopt)
	{
		pvt.tag = this;
		rosenbrock.findoptimum(pvt, pointopt, fvalopt);
	}


	//запись координат атомов класетра в файл
	void write_tetstofile(string filename)
	{
		writetetra3(filename, tets, ntets);
	}

	
	int copytoglobaltetra()
	{
		for(int i = 0; i < ntets; ++i)
			tetra[i] = tets[i];

		return ntets;
	}
};



//функция расчета R-фактора для его минимизации
float func_forrosenbrock(const nVector &point, void *tag)
{
	BuildClustConf &bc = *(static_cast<BuildClustConf*>(tag));
	
	bc.build_cluster_wpvtvector(point);

	//bc.write_tetstofile("temp.xyz");

	float rff = bc.dxrcl_cluster.rfactor_wcalcdxrcl(bc.tets, bc.ntets);

	return rff; //bc.dxrcl_cluster.rfactor_wcalcdxrcl(bc.tets, bc.ntets);
}


