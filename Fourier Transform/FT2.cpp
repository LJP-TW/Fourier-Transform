#include "FT.h"

using namespace std;

void FT::LowpassFilter(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;
	double**pFreq = new double*[M];
	for (int newcnt = 0;newcnt < M;newcnt++)
	{
		pFreq[newcnt] = new double[N];// �ť߸��W�v�}�C
	}
	for (int forzero_i = 0;forzero_i < M;forzero_i++)//��l��
	{
		for (int forzero_j = 0;forzero_j < N;forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	double cutoff = w / 16;
	double n = 4.0;
	for (int i = 0;i < M;i++)
	{
		for (int j = 0;j < N;j++)
		{
			double Huv = 1.0 / (1.0 + pow(hypot<double, double>(j-h/2, i-w/2) / cutoff, 2.0*n));
			//pFreq[i][j] *= Huv;
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], 2)*Huv + pow(FreqImag[i][j], 2)*Huv);
			OutputImage[i][j] = pFreq[i][j];
		}
	}
	for (int delcnt = 0;delcnt < M;delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[]pFreq;
    //cout << "LowpassFilter" << endl;
}

void FT::HighpassFilter(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;
	double**pFreq = new double*[M];
	for (int newcnt = 0;newcnt < M;newcnt++)
	{
		pFreq[newcnt] = new double[N];// �ť߸��W�v�}�C
	}
	for (int forzero_i = 0;forzero_i < M;forzero_i++)//��l��
	{
		for (int forzero_j = 0;forzero_j < N;forzero_j++)
		{
			pFreq[forzero_i][forzero_j] = 0.0f;
		}
	}
	//-------------------------------------------
	double cutoff = w / 8;
	double n = 4.0;
	for (int i = 0;i < M;i++)
	{
		for (int j = 0;j < N;j++)
		{
			double Huv = 1.0-(1.0 / (1.0 + pow(hypot<double, double>(j-h/2, i-w/2) / cutoff, 2.0*n)));
			//pFreq[i][j] *= Huv;
			pFreq[i][j] = sqrt(pow(FreqReal[i][j], 2)*Huv + pow(FreqImag[i][j], 2)*Huv);
			OutputImage[i][j] = pFreq[i][j];
		}
	}
	for (int delcnt = 0;delcnt < M;delcnt++)
	{
		delete[] pFreq[delcnt];
	}
	delete[]pFreq;
    cout << "HighpassFilter" << endl;
}
