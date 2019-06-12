#include "MyForm.h"
using namespace System;
using namespace System::Windows::Forms;
using namespace FourierTransform;

System::Void MyForm::lowpassFilterToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
	
	int w = dataManager->GetImageWidth();
	int h = dataManager->GetImageHeight();

	fourierTransformMethod->LowpassFilter(dataManager->GetInputImage(), dataManager->GetOutputImage(), dataManager->GetFreqReal(), dataManager->GetFreqImag(), h, w);
	Bitmap^ LPFImage = gcnew Bitmap(w, h);
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			int valuePixeli = dataManager->GetOutputImage()[i][j];
			if (valuePixeli > 255)
			{
				valuePixeli = 255;
			}
			else if (valuePixeli < 0)
			{
				valuePixeli = 0;
			}
			LPFImage->SetPixel(j, i, Color::FromArgb(valuePixeli, valuePixeli, valuePixeli));
		}
	}
	pictureBox_OutputImage->Image =LPFImage;
}

System::Void MyForm::highpassFilterToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
	int w = dataManager->GetImageWidth();
	int h = dataManager->GetImageHeight();

	fourierTransformMethod->HighpassFilter(dataManager->GetInputImage(), dataManager->GetOutputImage(), dataManager->GetFreqReal(), dataManager->GetFreqImag(), h, w);
	Bitmap^ HPFImage = gcnew Bitmap(w, h);
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			int valuePixeli = dataManager->GetOutputImage()[i][j];
			if (valuePixeli > 255)
			{
				valuePixeli = 255;
			}
			else if (valuePixeli < 0)
			{
				valuePixeli = 0;
			}
			HPFImage->SetPixel(j, i, Color::FromArgb(valuePixeli, valuePixeli, valuePixeli));
		}
	}
	pictureBox_OutputImage->Image = HPFImage;
}
