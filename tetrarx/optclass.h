#pragma once

#include "nvector.h"

class OPTIMIZE
{
public:
	//��������. countiterations - ��������� ���-�� ��������� ��������, countsucctrials - ���-�� �������� �������
	int countiterations, countsucctrials; 
	
	//�������� ��������� ������
	//stopbyepsfunc - �������� ����� ����� � ���������� ��������� ������� ������ epsfunc
	//stopbyepssteps - ��� ����� ���� ������ ������ epssteps
	int stopbyepsfunc, stopbyepssteps, stopbymaxiteration; 
	float fopt, x1, x2, x3;

	void *tag;

	OPTIMIZE(){n=0;};
	//�����������
	OPTIMIZE(
		int n_, 
		float epsfunc_, float epssteps_, 
		const nVector &steps_, float k_alpha_, float k_beta_,
		const nVector &pointmin_, const nVector &pointmax_, int isusepenalty_,
		const int maxitertation_,
		float (*func_)(const nVector &p, void *tag)
		);


	void setparameters(
	int n_, 
	float epsfunc_, float epssteps_, 
	const nVector &steps_, float k_alpha_, float k_beta_,
	const nVector pointmin_, const nVector pointmax_, int isusepenalty_,
	const int maxitertation_,
	float (*func_)(const nVector &p, void *tag)
	);

	//�-�� ��������� ���������� ������
	void     setalphabeta(float k_alpha_, float k_beta_	  ) { k_alpha  = k_alpha_;  k_beta   = k_beta_   ; };
	void      setepsilons(float epsfunc_, float epssteps_ ) { epsfunc  = epsfunc_;  epssteps = epssteps_ ; };
	void  setpointbarrier(const nVector &pointmin_, 
		                  const nVector &pointmax_        ) { pointmin = pointmin_; pointmax = pointmax_ ; };
	void  setmaxiteration(int maxiteration_               ) { maxiteration = maxiteration_               ; };
	void         setsteps(const nVector &steps0_          ) { steps0 = steps0_; steps = steps0_          ; };
	void  setminimizefunc(float (*func_)(const nVector &p,
						  void *tag)				  ) { func = func_								 ; };
	void  setisusepenalty(int isusepenalty_               ) { isusepenalty = isusepenalty_		         ; };
	void setiterationmode(int iterationmode_              ) { iterationmode = iterationmode_	         ; };

	//��������������� �-�� ������ �������� pointopt, fvalopt - ��������� ����� � �������� ��������
	void findoptimum(const nVector &point0, nVector &pointopt, float &fvalopt);

private:
	//n - ���������� ����������(����������, ���������)
	int     n;
	
	//��� �������� 0-����� �� epsfunc, epssteps �/��� maxiteration, 1-����� �� ������ ��������� �������� maxiteration
	//isusepenalty - ������������ �������� ������� - 1, �� ������������ - 0
	int     isusepenalty, maxiteration, iterationmode;  
	
	//epsilon - �������� �.�. ���������, ������ ��� ������� ����������� ����� � ���������� ������
	float   epsfunc, epssteps; 
	
	//k_alpha - �����.(>1 ����������=3) ����������(��������� �� ��� � �������� �����������) 
	//k_beta  - �����. ������ (�������������, ���. �� ��� � ���. �����������)    0 < abs(k_beta) < 1
	float   k_alpha, k_beta; 
	
	//���. ���/��������� ���/�������� ���������� ���������� �����/���.����������
	nVector steps0, steps, pointmin, pointmax, bvc;
	
	//������ �� �-�� ��� ���. ������ �������
	float   (*func)(const nVector &p, void *tag);

	//���. ������ ��������� ������ ����� � i-�����������
	int *attr;
	
	//�������, e-����������������� �����
	nVector *a, *e, *d;
	
	//�������� ����������, �� ���������� � �������� ��������� pointmin-pointmax
	int  ispointinborders(const nVector &p);

	//������� �� ������� ��� ������ ���������� �� ������� ��������� pointmin-pointmax
	float func_wpenalty(const nVector &p, float fvalold);
	
	//�������� ��������. ��������. ������ - (1,0,0) (0,1,0) (0,0,1) �������� ��� 3� ����������
	void setbazis();

	//��������� ��������������� �����-������
	void grammshmithortoganalization();

	//����� � �������. ������������ ������� ����������
	void rosenbrockprocedure(const nVector &point0, nVector &point_new,  float fval0, float &fval_new);
	
public:
	~OPTIMIZE()
	{
		delete [] a;
		delete [] e;
		delete [] d;
		delete [] attr;
	}

};





OPTIMIZE::OPTIMIZE(
	int n_, 
	float epsfunc_, float epssteps_, 
	const nVector &steps_, float k_alpha_, float k_beta_,
	const nVector &pointmin_, const nVector &pointmax_, int isusepenalty_,
	const int maxitertation_,
	float (*func_)(const nVector &p, void *tag))
{
	setparameters( n_, epsfunc_, epssteps_, steps_, k_alpha_, k_beta_, 
		pointmin_, pointmax_, isusepenalty_, maxitertation_, func_);
}



void OPTIMIZE::setparameters(
	int n_, 
	float epsfunc_, float epssteps_, 
	const nVector &steps_, float k_alpha_, float k_beta_,
	const nVector pointmin_, const nVector pointmax_, int isusepenalty_,
	const int maxitertation_,
	float (*func_)(const nVector &p, void *tag)
	)
{
	if(n > 0)
	{
		delete [] a;
		delete [] e;
		delete [] d;
		delete [] attr;
	}


	n = n_;
	steps0.initnew(n);
	steps.initnew(n);
	pointmin.initnew(n);
	pointmax.initnew(n);
	bvc.initnew(n);
	
	a = new nVector[n];  //������� �������� ����� ���
	e = new nVector[n];  //����
	d = new nVector[n];  //���� ���������� �� ���

	for(int i = 0; i < n; ++i)
	{
	    a[i].initnew(n);
		e[i].initnew(n);
		d[i].initnew(n);
	}

	attr  = new int[n];  //�������� �������� �����
	memset(attr, 0, n * sizeof(int));

	setepsilons(epsfunc_, epssteps_);
	setalphabeta(k_alpha_, k_beta_);
	setsteps(steps_);
	setpointbarrier(pointmin_, pointmax_);
	setmaxiteration(maxitertation_);
	setminimizefunc(func_);
	setisusepenalty(isusepenalty_);
	setbazis();

	iterationmode = 0;
	countiterations = 0;
	countsucctrials = 0;
	stopbyepsfunc   = 0;
	stopbyepssteps  = 0;
	stopbymaxiteration = 0;


}


//��������� ���. ������
void OPTIMIZE::setbazis()
{
	//�������������� ����� 0...1
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			if(i == j) e[i].x[j] = 1.0f;
		}
	}
}

//�������� ��������� ����� �� ���������� � ��������� pointmin-pointmax
int OPTIMIZE::ispointinborders(const nVector &p)
{
	for(int i = 0; i < n; ++i)
		if(p.x[i] < pointmin.x[i] || p.x[i] > pointmax.x[i]) return -1;

	return 0;
}

//�-�� �� �������
float OPTIMIZE::func_wpenalty(const nVector &p, float fvalold)
{
	if(-1 == ispointinborders(p))
		return fvalold + 100.0f;
	else
		return (*func)(p, tag);
}

//��������� ��������������� �����-������
void OPTIMIZE::grammshmithortoganalization()
{
	//�������� ����� e - ��������. �����

	e[0] = a[0] * (1.0f / a[0].norm());

	for(int i = 1; i < n; ++i)
	{
		bvc = 0.0f;
		
		for(int k = 0; k < i; ++k) 
			bvc = bvc + e[k] * (a[i] * e[k]);

		bvc = a[i] - bvc;
		e[i] = bvc * (1.0f / bvc.norm());
	}
}

//��������� ������ ����������
void OPTIMIZE::rosenbrockprocedure(const nVector &point0, nVector &point_new,  float fval0, float &fval_new)
{
	memset(attr, 0, n * sizeof(int));
	point_new = point0;
	float fval_prev = fval0;

	for(int i = 0; i < n; ++i)
	{
		a[i] = 0.0f;
		d[i] = 0.0f;

		if(isusepenalty == 1)
			fval_new =  func_wpenalty(point_new + e[i] * steps.x[i], fval_prev);
		else
			fval_new = (*func)(point_new + e[i] * steps.x[i], tag);

		if(fval_new >= fval_prev){         //�������
			attr[i] = 1;
			steps.x[i] = steps.x[i] * k_beta;
		} 
		else{
			d[i] =  e[i] * steps.x[i];
			point_new = point_new + d[i];
			steps.x[i] = steps.x[i] * k_alpha;
			fval_prev = fval_new;
			++countsucctrials;
		}
	}

	//���������� �������� �� ������� ����� �������� �����
	//����������������� ����� � ��������� ����� � ����� ������������
	for(int i = 0; i < n; ++i)
		if(attr[i] > 0)	a[i] = e[i];
		else for(int j = i; j < n; ++j) a[i] = a[i] + d[j];
}

//��������������� ����� �������� � ���������
void OPTIMIZE::findoptimum(const nVector &point0, nVector &pointopt, float &fvalopt)
{
	countiterations = 0;
	countsucctrials = 0;
	stopbyepsfunc   = 0;
	stopbyepssteps  = 0;
	stopbymaxiteration = 0;;
	
	tag = point0.tag;

	steps = steps0;  //��������. ����
	setbazis();      //������� ����� 0..1

	nVector point_prev(n), point_new(n);
	point_prev = point0;

	float fval_prev;
	float fval_new;

	fval_prev = (*func)(point_prev, tag);  //���. �-�� � ��������� �����

	if(isusepenalty == 1)
		fval_prev =  func_wpenalty(point_prev, fval_prev);

	for(int i = 0; i < maxiteration; ++i)
	{
		rosenbrockprocedure(point_prev, point_new, fval_prev, fval_new); 

		//x1 = point_new.x[0];
		//x2 = point_new.x[1];
		//float ddd = fabs(fval_new - fval_prev);

		if(iterationmode == 0)  //���� ����� ������ �� epsilon
		{
			if(fabs(fval_new - fval_prev) < epsfunc) { stopbyepsfunc = 1;  break; }
			int countstepsdelta = 0;
			for(int i = 0; i < n; ++i) if(fabs(steps.x[i]) < epssteps) ++countstepsdelta;
			if(countstepsdelta == n){ stopbyepssteps = 1;  break;}
		}

		if(fval_new < fval_prev)  //���� ���.��. ������ ����.(�����)
		{
			point_prev = point_new;
			fval_prev = fval_new;
		}

		//��������������� ������-������
		grammshmithortoganalization();

		++countiterations;

	}

	if(countiterations == maxiteration)
		stopbymaxiteration = 1;

	pointopt = point_prev;
	pointopt.tag = point_prev.tag;
	fvalopt  = fval_prev;

	fopt = fvalopt;
	x1 = point_prev.x[0];
	x2 = point_prev.x[1];
	
	if(n == 3)
		x3 = point_prev.x[2];



}


