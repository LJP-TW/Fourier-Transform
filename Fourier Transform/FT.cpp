#include "FT.h"
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

FT::FT()
{
}

void FT::DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			DFT(FreqReal, FreqImag, InputImage,M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			OutputImage[i][j] = sqrt(pow(FreqReal[i][j], (double) 2.0) + pow(FreqImag[i][j], (double) 2.0));
		}
	}
	//-------------------------------------------
}

void FT::DFT(double ** pFreqReal, double ** pFreqImag, int ** InputImage, int h, int w, int u, int v)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;

	for (int y = 0; y < M; y++)
	{
		for (int x = 0; x < N; x++)
		{
			// 可先計算Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleDFT = (-1.0f * 2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleDFT);
			double s = sin(angleDFT);

			// 利用Eular's equation計算傅立葉之實虛數部分
			pFreqReal[u][v] += (double)InputImage[y][x] * c;
			pFreqImag[u][v] -= (double)InputImage[y][x] * s;
		}
	}

	pFreqReal[u][v] = pFreqReal[u][v] / (double)(M);
	pFreqImag[u][v] = pFreqImag[u][v] / (double)(M);
}

void FT::InverseDiscreteFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w)
{
	int M = h;
	int N = w;

	double** InverseReal = new double*[M];
	double** InverseImag = new double*[M];
	double** pFreq = new double*[M];

	for (int i = 0; i<M; i++)
	{
		InverseReal[i] = new double[N];
		InverseImag[i] = new double[N];
		pFreq[i] = new double[N]; // 傅立葉頻率陣列
	}

	for (int i = 0; i<M; i++)
	{
		for (int j = 0; j<N; j++)
		{
			InverseReal[i][j] = 0.0f;
			InverseImag[i][j] = 0.0f;
			pFreq[i][j] = 0.0f;
		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			InverseDFT(InverseReal, InverseImag,FreqReal, FreqImag, M, N, j, i);
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			// 將計算好的傅立葉實數與虛數部分作結合 
			pFreq[i][j] = sqrt(pow(InverseReal[i][j], (double) 2.0) + pow(InverseImag[i][j], (double) 2.0));
			// 結合後之頻率域丟入影像陣列中顯示 
			OutputImage[i][j] = pFreq[i][j];
			//存下反傅立葉實數與虛數部分
			FreqReal[i][j] = InverseReal[i][j];
			FreqImag[i][j] = InverseImag[i][j];

		}
	}
	//-------------------------------------------
	for (int i = 0; i < M; i++)
	{
		delete[] pFreq[i];
		delete[] InverseReal[i];
		delete[] InverseImag[i];

	}
	delete[] pFreq;
	delete[] InverseReal;
	delete[] InverseImag;

}

void FT::InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y)
{
	// M = N 必須是方陣
	int M = h;
	int N = w;

	for (int v = 0; v < M; v++)
	{
		for (int u = 0; u < N; u++)
		{
			// 可先計算Eular's equation e^{j*theta} = cos{theta}+j*sin{theta}			
			double angleIDFT = (2.0f * 3.14159 * (double)(u*x + v*y) / (double)M);
			double c = cos(angleIDFT);
			double s = sin(angleIDFT);

			// 利用Eular's equation計算傅立葉之實虛數部分
			InverseReal[x][y] += (pFreqReal[v][u] * c - pFreqImag[v][u] * s);
			InverseImag[x][y] += (pFreqReal[v][u] * s + pFreqImag[v][u] * c);
		}
	}
	InverseReal[x][y] = InverseReal[x][y] / (float)M;
	InverseImag[x][y] = InverseImag[x][y] / (float)M;
}

void FT::FastFourierTransform(int** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int N)
{
#ifdef _DEBUG
    cout << "FastFourierTransform" << endl;
#endif

    //-------------------------------------------
    // 做成 complex
    vector<vector<complex<double>>> x(N);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            complex<double> c(InputImage[i][j], 0);
            x[i].push_back(c);
        }
    }

    // 做 x 軸方向的 FFT
    for (int l = 0; l < N; ++l)
    {
        /* bit-reversal permutation */
        for (int i = 1, j = 0; i < N; ++i)
        {
            for (int k = N >> 1; !((j ^= k)&k); k >>= 1);
            if (i > j) swap(x[l][i], x[l][j]);
        }

        /* dynamic programming */
        for (int k = 2; k <= N; k <<= 1)
        {
            float t = -2.0 * M_PI / k;
            complex<double> dw(cos(t), sin(t));

            // 每k個做一次FFT
            for (int j = 0; j < N; j += k)
            {
                // 前k/2個與後k/2的三角函數值恰好對稱，
                // 因此兩兩對稱的一起做。
                complex<double> w(1, 0);
                for (int i = j; i < j + k / 2; i++)
                {
                    complex<double> a = x[l][i];
                    complex<double> b = x[l][i + k / 2] * w;
                    x[l][i] = a + b;
                    x[l][i + k / 2] = a - b;
                    w *= dw;
                }
            }
        }
    }

    // 做 y 軸方向 FFT
    for (int l = 0; l < N; ++l)
    {
        /* bit-reversal permutation */
        for (int i = 1, j = 0; i < N; ++i)
        {
            for (int k = N >> 1; !((j ^= k)&k); k >>= 1);
            if (i > j) swap(x[i][l], x[j][l]);
        }

        /* dynamic programming */
        for (int k = 2; k <= N; k <<= 1)
        {
            float t = -2.0 * M_PI / k;
            complex<double> dw(cos(t), sin(t));

            // 每k個做一次FFT
            for (int j = 0; j < N; j += k)
            {
                // 前k/2個與後k/2的三角函數值恰好對稱，
                // 因此兩兩對稱的一起做。
                complex<double> w(1, 0);
                for (int i = j; i < j + k / 2; i++)
                {
                    complex<double> a = x[i][l];
                    complex<double> b = x[i + k / 2][l] * w;
                    x[i][l] = a + b;
                    x[i + k / 2][l] = a - b;
                    w *= dw;
                }
            }
        }
    }
    
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            FreqReal[i][j] = x[i][j].real() / N;
            FreqImag[i][j] = x[i][j].imag() / N;

            // 將計算好的傅立葉實數與虛數部分作結合 
            OutputImage[i][j] = sqrt(pow(FreqReal[i][j], (double) 2.0) + pow(FreqImag[i][j], (double) 2.0));
        }
    }
    //-------------------------------------------
}

void FT::InverseFastFourierTransform(int ** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int N)
{
#ifdef _DEBUG
    cout << "InverseFastFourierTransform" << endl;
#endif
    //-------------------------------------------
    // 做成 complex
    vector<vector<complex<double>>> x(N);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            complex<double> c(FreqReal[j][i], FreqImag[j][i]);
            x[i].push_back(c);
        }
    }

    // 做 x 軸方向的 FFT
    for (int l = 0; l < N; ++l)
    {
        /* bit-reversal permutation */
        for (int i = 1, j = 0; i < N; ++i)
        {
            for (int k = N >> 1; !((j ^= k)&k); k >>= 1);
            if (i > j) swap(x[l][i], x[l][j]);
        }

        /* dynamic programming */
        for (int k = 2; k <= N; k <<= 1)
        {
            float t = 2.0 * M_PI / k;
            complex<double> dw(cos(t), sin(t));

            // 每k個做一次FFT
            for (int j = 0; j < N; j += k)
            {
                // 前k/2個與後k/2的三角函數值恰好對稱，
                // 因此兩兩對稱的一起做。
                complex<double> w(1, 0);
                for (int i = j; i < j + k / 2; i++)
                {
                    complex<double> a = x[l][i];
                    complex<double> b = x[l][i + k / 2] * w;
                    x[l][i] = a + b;
                    x[l][i + k / 2] = a - b;
                    w *= dw;
                }
            }
        }
    }

    // 做 y 軸方向 FFT
    for (int l = 0; l < N; ++l)
    {
        /* bit-reversal permutation */
        for (int i = 1, j = 0; i < N; ++i)
        {
            for (int k = N >> 1; !((j ^= k)&k); k >>= 1);
            if (i > j) swap(x[i][l], x[j][l]);
        }

        /* dynamic programming */
        for (int k = 2; k <= N; k <<= 1)
        {
            float t = 2.0 * M_PI / k;
            complex<double> dw(cos(t), sin(t));

            // 每k個做一次FFT
            for (int j = 0; j < N; j += k)
            {
                // 前k/2個與後k/2的三角函數值恰好對稱，
                // 因此兩兩對稱的一起做。
                complex<double> w(1, 0);
                for (int i = j; i < j + k / 2; i++)
                {
                    complex<double> a = x[i][l];
                    complex<double> b = x[i + k / 2][l] * w;
                    x[i][l] = a + b;
                    x[i + k / 2][l] = a - b;
                    w *= dw;
                }
            }
        }
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            FreqReal[i][j] = x[i][j].real() / N;
            FreqImag[i][j] = x[i][j].imag() / N;

            // 將計算好的傅立葉實數與虛數部分作結合 
            OutputImage[i][j] = sqrt(pow(FreqReal[i][j], (double) 2.0) + pow(FreqImag[i][j], (double) 2.0));
        }
    }
    //-------------------------------------------
}
