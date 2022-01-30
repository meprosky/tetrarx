#pragma once

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<algorithm>
#include<cfloat>

using namespace std;


class DXRCL
{
	
	
public:
	float *s_exp, *hs_exp, *hs_exp_normby1, *hs_calc, *hs_calc_normby1;
	
	float *Io;
    float *Inorm;
	float *Hs;

	float rfactor;
	float sum_hs_exp, sum_hs_exp_normby1;			//сумма для вычисления R-фактора, вычислется в функции, где читаются данные

	
private:
	string comments1, comments2,	    //две строки заголовка вх. файла
		comments3, comments4;			//две строки для вых. файла

	int Na, Ns;							//число атомов, число сортов
	float lambdaXray;					//длина волны рентг. излучения

	int Nparc;							//число парциальных распределений пар атомов, выводимых 
	//в файлы для "ху". Если nparcr = 0, то не выводим и след. строка
	//массив parc должен отсутствовать

	int *parc;							//массив пар атомов

	float 
		Smin, Smax, dS, 
		Alpha, Drmax, delR, R0e, eps;  //Smin, Smax, dS, Al, Rmax, dR, R0e, eps

	float CHFE;							//число формульных единиц в кластере

	int   *ps;							//массив сортов атомов

	float *x, *y, *z;					//координаты атомов x, y, z

	int Nir;							//число интервалов по R
	int NSma;							//число интервалов по S
	float *Rij_0, *DispRij_0;			//"сырые" данные сумм R, квадартов
	int *NRij_0;						//кол-ва
	float *Rij, *DispRij;				//средние R, дисперсии
	int *NRij;							//кол-ва по R
	int *SA1, *SA2;						//вспомог. массивы индексов
	int NRma;							//кол-во интервалов для D(r)
	int Nreps, Npairs;					//кол-вл интервалов по eps(расстояний), пар перестановок

	float *DF1, *DF2;					 //первая и вторая дисперсионные поправки

	//9 коэффициентов на каждый сорт атома 
	//для расчета функций атомного рассеяния 
	float 
		*FA1, *FA2, *FA3, *FA4, 
		*FB1, *FB2, *FB3, *FB4, *FC, *FE;


	float *Gf2_expSAlpha;           //квадраты ф-ий обостряющего G-фактора в зав-ти от S память выделяется в функции где расчитывается
	float *KFAR;					//коэф-ты ФАР  в зав-ти от S,  память выделяется в функции где расчитывается
	float GfactorSum1;

	//************************************************************************
	//     Вспомогательная матрица индексов для доступа к элементам массивов
	//************************************************************************
	int Ns1;
	int *idxA;
	int *invidxSA1;
	int *invidxSA2;

	int Nall;  // Nall = Npairs * Nreps;          //длина массивов с учетом идентичности пар (1-2 = 2-1)
	                                              //например для 3 сортов атомов будет 1-1 1-2 1-3 2-2 2-3 3-3 всего 6 пар

	
public:
		
	DXRCL()
	{
	}

	DXRCL(string filename)
	{
		dxrclinit(filename);
	}

	//Конструктор
	void dxrclinit(string filename)
	{
		readdatadxrcl(filename);

		KFAR  = new float[NSma * Ns * Ns];      //массив ФАР - суть 3х мерный массив, строки - значения S 
		//и соотв. ему 2х мерный массив ФАР i и j атома
		//возвращаемое значение

		Gf2_expSAlpha  = new float[NSma];  //возвращаемое значение

		calcFAR();

		Io    = new float[NSma];
		Inorm = new float[NSma];
		Hs    = new float[NSma];

		Npairs = Ns * (Ns + 1) / 2;				// число перестановок сортов пар атомов

		//************************************************************************
		//     Вспомогательная матрица индексов для доступа к элементам массивов
		//************************************************************************
		Ns1 = Ns + 1;
		idxA      = new int[Ns1*Ns1];
		invidxSA1 = new int[Npairs];
		invidxSA2 = new int[Npairs];

		int cc = 0;
		for(int i = 1; i <= Ns; ++i)         //заполняем матрицу индексов для доступа к элементам массивов
		{
			for(int j = i; j <=Ns; ++j)
			{
				idxA[i*Ns1+j] = cc; 
				idxA[j*Ns1+i] = cc;
				invidxSA1[cc] = i;
				invidxSA2[cc] = j;

				++cc;
			}
		}

	}

	~DXRCL()
	{
		delete [] KFAR;
		delete [] Gf2_expSAlpha;
		delete [] s_exp;    
	    delete [] hs_exp;
		delete [] hs_exp_normby1;
		delete [] hs_calc;
		delete [] hs_calc_normby1;

		delete [] Io;
		delete [] Inorm;
		delete [] Hs;
		//delete [] Dr;

		delete [] idxA;

		delete [] invidxSA1;
		delete [] invidxSA2;

		delete [] DF1;  
		delete [] DF2;  
		delete [] FA1;
		delete [] FA2;
		delete [] FA3;
		delete [] FA4;
		delete [] FB1;
		delete [] FB2;
		delete [] FB3;
		delete [] FB4;
		delete [] FC;
		delete [] FE;

	}

	
private:

	//************************************************************************
	//
	//     Чтение входных данных
	//
	//************************************************************************
	inline void readdatadxrcl(string filename)
	{
		ifstream inputDataFile(filename);
		stringstream strstream;


		if(!inputDataFile)
		{
			cout << "Ошибка чтения входного файла" << endl;
			getchar();
			exit(1);
		}

		strstream << inputDataFile.rdbuf();  //весь файл в строковый поток
		inputDataFile.close();               //вх. файл больше не нужен, закрываем 

		string inputdata = strstream.str();  //весь поток направляем в строку

		getline(strstream, comments1); //читаем строку 1
		getline(strstream, comments2); //читаем строку 2

		strstream >> Na >> Ns; //число атомов, число сортов

		strstream >> lambdaXray; // длина волны рентг. излучения

		DF1 = new float[Ns];  //первая дисперсионная поправка
		DF2 = new float[Ns];  //вторая дисперсионная поправка

		for(int i = 0; i < Ns; ++i) strstream >> DF1[i];
		for(int i = 0; i < Ns; ++i) strstream >> DF2[i];

		FA1 = new float[Ns];  //9 коэффициентов на каждый сорт атома 
		FA2 = new float[Ns];  //для расчета функций атомного рассеяния 
		FA3 = new float[Ns];
		FA4 = new float[Ns];
		FB1 = new float[Ns];
		FB2 = new float[Ns];
		FB3 = new float[Ns];
		FB4 = new float[Ns];
		FC  = new float[Ns];

		for(int i = 0; i < Ns; ++i) strstream >> FA1[i];
		for(int i = 0; i < Ns; ++i) strstream >> FA2[i];
		for(int i = 0; i < Ns; ++i) strstream >> FA3[i];
		for(int i = 0; i < Ns; ++i) strstream >> FA4[i];
		for(int i = 0; i < Ns; ++i) strstream >> FB1[i];
		for(int i = 0; i < Ns; ++i) strstream >> FB2[i];
		for(int i = 0; i < Ns; ++i) strstream >> FB3[i];
		for(int i = 0; i < Ns; ++i) strstream >> FB4[i];
		for(int i = 0; i < Ns; ++i) strstream >> FC[i];

		FE = new float[Ns];  //"формульная единица" или единица состава 
		// Например, Al2O3 - FE(1) = 2., FE(2) = 3. 
		//( Признак сорта атома PS -  Al - 1, O - 2) 

		for(int i = 0; i < Ns; ++i) strstream >> FE[i];

		strstream >> Nparc; //число парциальных распределений пар атомов, выводимых 
		//в файлы для "ху". Если nparcr = 0, то не выводим и след. строка
		//должна отсутствовать


		if(Nparc > 0)
		{
			parc = new int[Nparc];
			string s;
			getline(strstream, s);
			getline(strstream, s);

			replace(s.begin(), s.end(), ',', ' ');

			stringstream ss2(s);

			for(int i = 0; i < Nparc; ++i) ss2 >> parc[i];

		}
		else if(Nparc < 0)
		{
			cout << "Ошибка числа парциальных распределений Nparc";
			exit(1);
		}
		else //0
		{
		}

		string s;
		getline(strstream, s);
		replace(s.begin(), s.end(), ',', ' ');

		stringstream ss2(s);
		ss2 >> Smin >> Smax >> dS >> Alpha >> Drmax >> delR >> R0e >> eps; //Smin, Smax, dS, Al, Rmax, dR, R0e, eps

		NSma = (int)((Smax - Smin) / dS + 1.1f);       //кол-во интервалов dS в диапазоне Smin-Smax (вычисляем здесь, позже понадобится)

		strstream >> CHFE; // число формульных единиц в кластере

		getline(strstream, s);
		getline(strstream, s);

		int NSma_frstr;

		strstream >> NSma_frstr;

		if(NSma != NSma_frstr)
		{
			cout << "Ошибка чтения интервалов по S";
			exit(1);
		}

		s_exp   = new float[NSma]; 
		hs_exp  = new float[NSma]; 
		hs_calc = new float[NSma]; 
		hs_exp_normby1  = new float[NSma]; 
		hs_calc_normby1 = new float[NSma];

		strstream >> s_exp[0] >> hs_exp[0];

		float hs_min = hs_exp[0], hs_max = hs_exp[0];

		for(int i = 1; i < NSma; ++i) 
		{
			strstream >> s_exp[i] >> hs_exp[i];

			float sss = fabs(s_exp[i] - s_exp[i-1]) - fabs(dS);
			float sss2 = s_exp[i];
			float sss3 = s_exp[i-1];
			float sss4 = sss2-sss3;
			float sss5 = sss4-dS;

			if(fabs(s_exp[i] - s_exp[i-1]) - fabs(dS) > 2 * FLT_EPSILON)
			{
				cout << "Ошибка. Проверьте dS в заголовке и графике H(s)";
				exit(1);
			}

			if(hs_exp[i] > hs_max) hs_max = hs_exp[i]; //ищем min-max чтобы нормировать
			if(hs_exp[i] < hs_min) hs_min = hs_exp[i];

		}

		if( fabs(s_exp[0] - Smin) > 2 * FLT_EPSILON || fabs(s_exp[NSma-1] - Smax) > 2 * FLT_EPSILON)
		{
			cout << "Ошибка. Проверьте Smin, Smax в заголовке и графике H(s)";
			exit(1);
		}

		float hs_norm = fabs(hs_max - hs_min);


		sum_hs_exp = 0.0f;
		sum_hs_exp_normby1 = 0.0f;


		//нормируем H(s) на макс.
		for(int i = 0; i < NSma; ++i) 
		{
			hs_exp_normby1[i] = hs_exp[i] / hs_norm;
			sum_hs_exp += fabs(hs_exp[i]);  //для вычисления R фактора в дальнейшем
			sum_hs_exp_normby1 += fabs(hs_exp_normby1[i]);
		}

		/*
		ps = new   int[Na];   //сорт атома
		x =  new float[Na];   //координаты атомов x, y, z
		y =  new float[Na];
		z =  new float[Na];

		for(int i = 0; i < Na; ++i) 	
		strstream >> ps[i] >> x[i] >> y[i] >> z[i];
		*/


	}

	inline float findRmax()
	{
		float minx = x[0], miny = y[0], minz = z[0];
		float maxx = x[0], maxy = y[0], maxz = z[0];

		for(int i = 0; i < Na; ++i) //ищем мин и макс значения координат
		{
			if(x[i] < minx) minx = x[i];
			if(y[i] < miny) miny = y[i];
			if(z[i] < minz)	minz = z[i];

			if(x[i] > maxx) maxx = x[i];
			if(y[i] > maxy) maxy = y[i];
			if(z[i] > maxz)	maxz = z[i];
		}

		float dx=maxx-minx, dy=maxy-miny, dz=maxz-minz;

		return sqrtf(dx*dx + dy*dy + dz*dz);            //макс. диагональ
	}

	//************************************************************************
	//     Вычисление совокупности межатомных расстояний в кластере
	//     с суммированием числа одинаковых (в пределах EPS) расстояний
	//     между одинаковыми парами атомов (A - B и B - A - одинаковы).
	//     Выч. среднее межат. расстояние и дисперсия этого расстояния
	//     по совокупности расстояний в пределах EPS.
	//************************************************************
	inline void interatomic()
	{
		int Nam1 = Na - 1;
		int Ns1  = Ns + 1;
		int idxa, idxeps, idxG;

		//хранение необходимо по столбцам, чтобы были в порядке возрастания R
		//в строках интервалы по eps,  в столбцах пары атомов

		float eps_inv = 1.0f / eps;
		float dr, dr_p2, dx, dy, dz;


		//****************************************************************************
		//
		//   Вычисление межатомных расстояний, квадратов для дисперсий,
		//   сумм расстоняий
		//
		//****************************************************************************
		for(int i = 0; i < Nam1; ++i)
		{
			for(int j = i + 1; j < Na; ++j)
			{

				dx = x[j] - x[i];
				dy = y[j] - y[i];
				dz = z[j] - z[i];

				dr_p2 = dx*dx + dy*dy + dz*dz;     //квадрат
				dr = sqrtf(dr_p2);

				idxa = idxA[ ps[i] * Ns1 + ps[j] ];

				idxeps = (int)(dr * eps_inv);    //индекc чтоб попасть в интервал n*eps - (n+1)*eps

				idxG = idxeps * Npairs + idxa;   //глобальный индекс в массивах по Eps

				Rij_0[idxG] += dr;               //сумма расстояний, потом вычислим среднее

				DispRij_0[idxG] += dr_p2;        //для выч. дисперсии			

				NRij_0[idxG] += 1;               //подсчет кол-ва интервалов расстояний

			}
		}

	}

	//****************************************************************************
	//
	//   Вычисление значений ФАР в зависимости от S и не зависимых от
	//   числа атомов в кластере, от распеределений по расстояниям, дисперсиям и т.п.
	//
	//****************************************************************************
	void calcFAR()         //возвращаемое значение
	{
		int Ns2 = Ns * Ns;

		float L4Pi = 1.0f / (4.0f * PI);
		float ST;

		GfactorSum1 = 0.0f;
		for(int i = 0; i < Ns; ++i)
		{
			GfactorSum1 = GfactorSum1 + FE[i] * (FC[i] + FA1[i] + FA2[i] + FA3[i] + FA4[i]);
		}
		float GfactorSum1_inv = 1 / GfactorSum1;           //получаем обратную величину 
		//потом будем умножать (это быстрее чем деление)


		float S = Smin;
		for(int k = 0; k < NSma; ++k)
		{
			ST = S * L4Pi;
			ST = ST * ST;

			int _idxs = k * Ns2;                   //всп. индекс в 3х мерном массиве

			float GfactorSum2 = 0.0f;

			for (int i = 0; i < Ns; ++i)
			{
				float KFARFi = 
					FC[i] +
					FA1[i] * exp(-FB1[i] * ST) +			
					FA2[i] * exp(-FB2[i] * ST) +
					FA3[i] * exp(-FB3[i] * ST) +
					FA4[i] * exp(-FB4[i] * ST) + 
					DF1[i];

				float KP1 = KFARFi * KFARFi + DF2[i] * DF2[i];

				int _idx = i * Ns;             //всп. индекс в 3х мерном массиве

				KFAR[_idxs + _idx + i] = KP1;  //сорт i-i

				//сумма для расчета обостряющего G фактора
				KP1 = sqrtf(KP1);
				GfactorSum2 = GfactorSum2 + FE[i] * KP1;

				for (int j = i+1; j < Ns; ++j)
				{
					float KFARFj = 
						FC[j] + 
						FA1[j] * exp(-FB1[j] * ST) +			
						FA2[j] * exp(-FB2[j] * ST) +
						FA3[j] * exp(-FB3[j] * ST) +
						FA4[j] * exp(-FB4[j] * ST) + 
						DF1[j];

					float KP2 = 
						KFARFi * KFARFj + DF2[i] * DF2[j];

					KFAR[_idxs + _idx + j] = KP2;           //элемент матрицы - сорт i-j
					KFAR[_idxs + j*Ns + i] = KP2;           //элемент матрицы - сорт j-i (в принципе можно убрать)

				}
			}

			float Gf = GfactorSum2 * GfactorSum1_inv;
			float sAlpha = S * Alpha;
			Gf2_expSAlpha[k] = expf(- sAlpha * sAlpha) / (Gf * Gf);		// эксп деленная на квадрат обостряющего G фактор 
			//предвычисляем для дальнейшего расчета H(s)

			//Hs[k] =            интерференционная функция
			//  S * (Inorm[k] - IndScat / CHFE) * Gf2_expSAlpha[k]

			S = S + dS;
		}



	}

	void calcNrAvrDisperse()
	{
		int c = 0;
		float Rr;

		for(int i = 0; i < Nreps; ++i)							//пробегаем по интервалам eps
		{
			for(int j = 0; j < Npairs; ++j)						//пробегаем по комбинациям пар атомов 1-1 1-2 1-3 2-2 2-3 и т.д.
			{
				int k = i * Npairs + j;
				int Nr = NRij_0[k];

				if(Nr > 1)
				{
					NRij[c] = Nr;

					Rr = Rij[c] = Rij_0[k] / Nr;				//среднее расстояние в интервале n*Eps - (n+1)*Eps

					DispRij[c] =								//дисперсия
						sqrtf(
						fabsf(DispRij_0[k] - (float)Nr*Rr*Rr) /
						((float)Nr * (Nr - 1))
						); 

					SA1[c] = invidxSA1[j];						//сорт атома 1
					SA2[c] = invidxSA2[j];						//сорт атома 1 

					++c;

				}
				else if(Nr == 1)								//если расстояние в интервале n*Eps - (n+1)*Eps всего одно то дисперсия = 0
				{
					NRij[c] = 1;				
					Rij[c] = Rij_0[k];
					DispRij[c] = 0.0f;
					SA1[c] = invidxSA1[j];
					SA2[c] = invidxSA2[j];
					++c;
				}
			}
		}


	

	}

	//*****************************************************************************************
	//
	//   Блок вычислений функций I(s), H(s), D(r)
	//   
	//   Замечание! Поскольку сорта атомов начинаются с 1 (1,2,3...и т.д.) а массивы все с нуля,
	//   то при доступе к элементу массива зависящему от сорта атома, от сорта атома мы отнимаем 1
	//
	//*****************************************************************************************
	void calcINTS(		
		float &hs_min,
		float &hs_max
		)
	{

		float S = Smin;
		int Ns2 = Ns * Ns;

		float hs_calc_min = FLT_MAX;
		float hs_calc_max = - FLT_MIN;

		for(int k = 0; k < NSma; ++k)
		{
			int _idxs = k * Ns2;                                   //часть индекса для доступа в массив ФАР

			//*****************************************************************************************
			//     Расчет независимого рассеяния
			//		         Na
			//        I(s) = SUM fk(s)*fk(s)
			//               k=1
			//*****************************************************************************************

			float IndScat = 0.0f;									//слагаемое независимого рассеяния

			for(int i = 0; i < Na; ++i)
			{
				int _is = ps[i] - 1;								//см. замечание в заголовке
				int index = _idxs + _is*Ns + _is;					//индекс в массиве посчитанных коэффициентов ФАР
				//для сорта атома ps и значения S. См. замечание в заголовке

				//int index = _idxs + (ps[i] - 1) * (Ns + 1);         

				IndScat = IndScat + KFAR[index];                   //слагаемое независимого рассеяния
			}

			//*****************************************************************************************
			//        Расчет интерференционного слагаемого : 
			//                 Na-1 Na              *       *        
			//        I(s) =2* SUM SUM 1/2 *[fi(s)*fj(s) + fi(s)*fj(s)]*
			//                 i=1 j=i+1             
			//
			//        *SIN(s*rij)/(s*rij)*EXP(-0.5*(drij*s)**2) 
			//*****************************************************************************************

			float InterfScat = 0.0f;                                //интерференционное слагаемое

			for(int i = 0; i < Nir; ++i)
			{
				int _is = SA1[i] - 1;
				int _js = SA2[i] - 1;

				int index = _idxs + _is*Ns + _js;                 //индекс в массиве посчитанных коэффициентов ФАР
				//для сортов атомов SA1 и SA2 и значения S
				float d = S * DispRij[i];

				InterfScat = InterfScat +                          //вычисляем интерференционное слагаемое 
					KFAR[index] *
					NRij[i] * sinf(S * Rij[i]) /
					(S * Rij[i]) *
					expf(-0.5f * d * d );
			}

			float SumScat = IndScat + InterfScat * 2;              //суммарная интенсивность I(s)

			Io[k] = SumScat;                                       //I(s)

			Inorm[k] = SumScat / CHFE;                             //I(s) нормированная на единицу состава CHFE

			Hs[k] =                                                //H(s) интерференционная функция
				S * (Inorm[k] - IndScat / CHFE) * 
				Gf2_expSAlpha[k];                                  // Gf2_expSAlpha[k] =
			//    exp(- (S*Alpha)^2) / (GfactorSum2 / GfactorSum1)^2 ;

			if(Hs[k] > hs_calc_max) hs_calc_max = Hs[k];
			if(Hs[k] < hs_calc_min) hs_calc_min = Hs[k];

			S = S + dS;                                            //увеличиваем S на интервал dS
		}

		hs_min = hs_calc_min;
		hs_max = hs_calc_max;
		
	}

	//*****************************************************************************************
	//
	//      Расчет D(r)
	//
	//*****************************************************************************************
	
	void calcDr(const float *Hs, float* &Dr)   //возвращаемое значение
	{

		float FSL = 2.0f * PI * PI * R0e * GfactorSum1;
		NRma = (int)(Drmax / delR + 1.1f);
		Dr = new float[NRma];                 //D(r)
		float S;

		float R = 0;
		for(int i = 0; i < NRma; ++i)
		{
			S = Smin;
			float Sum = 0.0f;
			for(int k = 0; k < NSma; ++k)
			{
				Sum = Sum + Hs[k] * sinf(S * R);
				S = S + dS;
			}
			Dr[i] = FSL * R + Sum * dS;
			R = R + delR;
		}

	}


public:
	void Write_Rij_Nrij_Dispij(string fileName)
	{
		ofstream outf(fileName);
		outf << comments1 << endl << comments2 << endl << comments3 << endl;
		outf << setw(7)  << "i-j";
		outf << setw(7)  << "№пп"; 
		outf << setw(7)  << "N";
		outf << setw(15) << "Rij";
		outf << setw(15) << "DispRij" << endl;

		ostringstream str("");	


		for(int i = 1; i <= Ns; ++i)
			for(int j = i; j <= Ns; ++j)
			{
				int l = 1;
				for(int k = 0; k < Nir; ++k)
				{
					if(SA1[k] == i &&  SA2[k] == j)  
					{
						str.str("");
						str << i << "-" << j;
						outf << setw(7)  << str.str();
						outf << setw(7)  << l;
						outf << setw(7)  << NRij[k];
						outf << setw(15) << setprecision(5) << fixed << Rij[k];
						outf << setw(15) << setprecision(5) << fixed << DispRij[k];
						outf << resetiosflags(ios::scientific | ios::fixed);
						outf << endl;

						++l;
					}
				}
				outf << endl;
			}

	}

	void Write_Is(string fileName, const float *Io)
	{
		ofstream fo(fileName);

		fo  << comments1 << endl
			<< comments2 << endl
			<< comments3 << endl
			<< "I(S)" << endl
			<< NSma << " 30 30 1 0" << endl; 

		float S = Smin;
		for(int k = 0; k < NSma; ++k)
		{
			//S = Smin + k * dS;
			fo << setw(8)  << setprecision(4) << fixed << S;
			fo << setw(15) << setprecision(4) << fixed << Io[k] << endl;
			S = S + dS;
		}
	}

	void Write_Isnorm(string fileName, const float *Inorm)



	{
		ofstream fo(fileName);

		fo  << comments1 << endl
			<< comments2 << endl
			<< comments3 << endl
			<< "Inorm(S)" << endl
			<< NSma << " 3 3 1 0" << endl; 

		float S = Smin;
		for(int k = 0; k < NSma; ++k)
		{
			//S = Smin + k * dS;
			fo << setw(8)  << setprecision(4) << fixed << S;
			fo << setw(15) << setprecision(4) << fixed << Inorm[k] << endl;
			S = S + dS;
		}
	}

	void Write_Hs(string fileName, const float *Hs)
	{
		ofstream fo(fileName);

		fo  << comments1 << endl
			<< comments2 << endl
			<< comments3 << endl
			<< "H(S)" << endl
			<< NSma << " 1 1 1 0" << endl; 

		float S = Smin;
		for(int k = 0; k < NSma; ++k)
		{
			//S = Smin + k * dS;
			fo << setw(8)  << setprecision(4) << fixed << S;
			fo << setw(15) << setprecision(4) << fixed << Hs[k] << endl;
			S = S + dS;
		}
	}

	void Write_Dr(string fileName, const float *Dr)
	{
		ofstream fo(fileName);

		fo  << comments1 << endl
			<< comments2 << endl
			<< comments3 << endl
			<< "D(r)" << endl
			<< NRma << " 15 15 1 0" << endl; 

		float R = 0;
		for(int k = 0; k < NRma; ++k)
		{
			fo << setw(8)  << setprecision(3) << fixed << R;
			fo << setw(15) << setprecision(3) << fixed << Dr[k] << endl;

			R = R + delR;
		}
	}


public:

	void get_globaltetra()
	{
		Na = ATOMPERTETRA * Ntetra;

		ps = new   int[Na];   //сорт атома
		x =  new float[Na];   //координаты атомов x, y, z
		y =  new float[Na];
		z =  new float[Na];

		int countatom = 0;
		for(int i = 0; i < Ntetra; ++i)
		{
			for(int j = 0; j < tetra[i].na; ++j)
			{
				if(1 != tetra[i].a[j].deleted)
				{
					ps[countatom] = tetra[i].a[j].sort;
					x[countatom] = tetra[i].a[j].gc.x;
					y[countatom] = tetra[i].a[j].gc.y;
					z[countatom] = tetra[i].a[j].gc.z;

					++countatom;
				}

			}
		}

		Na = countatom;						//уточняем кол-во атомов
	}

	void get_configtetra(const Polyhedron *tets, int ntets)
	{
		Na = ATOMPERTETRA * ntets;

		ps = new   int[Na];   //сорт атома
		x =  new float[Na];   //координаты атомов x, y, z
		y =  new float[Na];
		z =  new float[Na];

		

		int countatom = 0;
		for(int i = 0; i < ntets; ++i)
		{
			for(int j = 0; j < tets[i].na; ++j)
			{
				if(1 != tets[i].a[j].deleted)
				{
					ps[countatom] = tets[i].a[j].sort;
					x[countatom] = tets[i].a[j].gc.x;
					y[countatom] = tets[i].a[j].gc.y;
					z[countatom] = tets[i].a[j].gc.z;

					if(countatom == 5 && (x[5] - x[0]) < 0.0000001)
					{
						cout << NULL;
					}

					++countatom;

					

					
				}

			}
		}

		Na = countatom;						//уточняем кол-во атомов
	}

	void dxrclcalc()
	{
		get_globaltetra();
		dxrclcalcmain();
	}

	void dxrclcalc(const Polyhedron *tets, int ntets)
	{
		get_configtetra(tets, ntets);
		dxrclcalcmain();
	}

	void dxrclcalcmain()
	{

		float Rmax = findRmax();			// Ищем макс. диагональ кластера

		Nreps = (int)(Rmax / eps + 1.1f);	//количество интервалов eps в макс. диагонали

		Nall = Npairs * Nreps;				//длина массивов с учетом идентичности пар (1-2 = 2-1)
											//например для 3 сортов атомов будет 1-1 1-2 1-3 2-2 2-3 3-3 всего 6 пар

		Rij_0   =    new float[Nall];     //межатомные расстояния
		DispRij_0  = new float[Nall];     //дисперсии расстояний
		NRij_0  =    new   int[Nall];     //кол-во расстояний в интервале eps

		unsigned int sz_float = Nall * sizeof(float);
		unsigned int sz_int   = Nall * sizeof(int);

		memset(Rij_0,      0, sz_float);   //обнуляем
		memset(DispRij_0,  0, sz_float);
		memset(NRij_0,     0, sz_int);


		//************************************************************************
		//		САМЫЕ ЗАТРАТНЫЕ ПО ВРЕМНИ ВЫЧИСЛЕНИЯ (ибо зависят от от N(N-1)/2 )
		//		Вычисление совокупности межатомных расстояний в кластере
		//		с суммированием числа одинаковых (в пределах EPS) расстояний
		//		между одинаковыми парами атомов (A - B и B - A - одинаковы).
		//		Выч. среднее межат. расстояние и дисперсия этого расстояния
		//		по совокупности расстояний в пределах EPS.
		//************************************************************************
		interatomic();


		//************************************************************************
		//	Подсчитываем число ненулевых элементов Nir
		//	т.е. тех межат. расстояния которых встречаются хотя бы один раз
		//	делается для выделения памяти для массивов точно под размер
		//************************************************************************
		Nir = 0; 
		for(int i = 0; i < Nreps; ++i){      //пробегаем по интервалам eps
			for(int j = 0; j < Npairs; ++j){ //пробегаем по комбинациям пар атомов 1-1 1-2 1-3 2-2 2-3 и т.д.
				if(NRij_0[i * Npairs + j]){  //не равно 0
					++Nir;                  //сколько ненулевых интервалов
				}
			}
		}

		//****************************************************************************
		//
		//   Вычисление средних расстояний, их количества в интервале, дисперсий 
		//   попутно "cжимаем" массив, т.е. исключаем элементы с NRij = 0 (нет расстояний)
		//
		//****************************************************************************
		Rij     = new float[Nir];     //межатомные расстояния
		DispRij = new float[Nir];     //дисперсии расстояний
		NRij    = new   int[Nir];     //кол-во расстояний в интервале eps
		SA1     = new   int[Nir];     //первый сорт атома в паре
		SA2     = new   int[Nir];     //второй сорт атома в паре
				
		calcNrAvrDisperse();
	
		delete [] Rij_0;                                       
		delete [] NRij_0;
		delete [] DispRij_0;

		//NSma = (int)((Smax - Smin) / dS + 1.1f);       //кол-во интервалов dS в диапазоне Smin-Smax (выч. еще при чтении данных из файла)

		//****************************************************************************
		//
		//	ПРЕДВЫЧИСЛЕНИЯ ФАР СДЕЛАНЫ В САМОМ НАЧАЛЕ РАБОТЫ ПРОГРАММЫ
		//	ЭТОТ БЛОК ЗДЕСЬ ПРОПУЩЕН
		//   Вычисление значений ФАР в зависимости от S и не зависимых от
		//   числа атомов в кластере, от распеределений по расстояниям, дисперсиям и т.п.
		//
		//****************************************************************************

		//*****************************************************************************************
		//
		//	 ОСНОВНОЙ БЛОК
		//   Вычисление функций I(s), H(s), D(r)
		//   
		//*****************************************************************************************
		
		float hs_min, hs_max;

		calcINTS(hs_min, hs_max);

		//нормируем H(s)
		float hs_norm = fabs(hs_max - hs_min);
		for(int i = 0; i < NSma; ++i) 
		{
			hs_calc[i]         = Hs[i];
			hs_calc_normby1[i] = Hs[i] / hs_norm;
		}

		//*****************************************************************************************
		//	ПОКА НЕ БУДЕМ ДЕЛАТЬ
		//      Расчет D(r)
		//*****************************************************************************************
		//float *Dr;
		//calcDr(Hs, Dr);

		comments3 = "tetrarx";

		//Write_Hs("HSTETRA.xy", Hs);
		//Write_Hs("HSTETRA_calc.xy", hs_calc);

		//освобождаем память
		delete [] Rij;                                      
		delete [] NRij;
		delete [] DispRij;
		delete [] SA1;                                      
		delete [] SA2;
		


	}

	float rfactorcalc()
	{
		float f = 0.0f;
		for(int i = 0; i < NSma; ++i) 
		{
			f += fabs(hs_exp_normby1[i] - hs_calc_normby1[i]);
		}
		rfactor = f / sum_hs_exp_normby1;
		return rfactor;
	}

	float rfactor_wcalcdxrcl()
	{
		dxrclcalc();
		return rfactorcalc();
	}

	float rfactor_wcalcdxrcl(const Polyhedron *tets, int ntets)
	{
		dxrclcalc(tets, ntets);
		return rfactorcalc();
	}



};