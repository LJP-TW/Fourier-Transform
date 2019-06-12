#pragma once
#include <iostream>
#include <vector>
#include <complex>
class FT
{
private:

public:
    FT();
    void DiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
    void DFT(double** pFreqReal, double** pFreqImag, int** InputImage, int h, int w, int u, int v);

    void InverseDiscreteFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
    void InverseDFT(double ** InverseReal, double ** InverseImag, double ** pFreqReal, double ** pFreqImag, int h, int w, int x, int y);

    void FastFourierTransform(int** InputImage, int ** OutputImage, double ** FreqReal, double ** FreqImag, int N);

    void InverseFastFourierTransform(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int N);

    void LowpassFilter(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);
    void HighpassFilter(int** InputImage, int** OutputImage, double ** FreqReal, double ** FreqImag, int h, int w);

private:

};



