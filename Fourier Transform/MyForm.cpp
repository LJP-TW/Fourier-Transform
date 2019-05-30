#include "MyForm.h"
using namespace System;
using namespace System::Windows::Forms;
using namespace FourierTransform;

System::Void MyForm::loadImageToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
    openFileDialog1->ShowDialog();
}

System::Void MyForm::openFileDialog1_FileOk(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e)
{
    //將影像讀入並設儲存至DataManager的資料結構中
    Bitmap^ OriginalImage = gcnew Bitmap(openFileDialog1->FileName);
    pictureBox_SourceImage->Image = OriginalImage;
    int scale = 1;
    int M = OriginalImage->Height * scale;
    int N = OriginalImage->Width * scale;

    dataManager = new DataManager(M, N);
    for (int i = 0; i < OriginalImage->Height; i++)
    {
        for (int j = 0; j < OriginalImage->Width; j++)
        {
            Color srcColor = OriginalImage->GetPixel(j, i); // 擷取每個點的顏色
            int srcGrey = srcColor.R*0.299 + srcColor.G*0.587 + srcColor.B*0.144; // 彩色三通道轉成灰階
            dataManager->SetPixel(j, i, srcGrey);
        }
    }

    std::cout << "-Image has been loaded-" << std::endl;
}

System::Void MyForm::discreteFourierTransformToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
    int w = dataManager->GetImageWidth();
    int h = dataManager->GetImageHeight();

    // 利用傅立葉之平移性，平移頻率
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            int valuePixeli = dataManager->GetInputImage()[i][j];
            valuePixeli = valuePixeli * pow((float)-1, (float)(i + j));
            dataManager->SetPixel(j, i, valuePixeli);
        }
    }

    //將算出頻率資訊傳入輸出影像
    fourierTransformMethod->DiscreteFourierTransform(dataManager->GetInputImage(), dataManager->GetOutputImage(), dataManager->GetFreqReal(), dataManager->GetFreqImag(), h, w);
    Bitmap^ DFTImage = gcnew Bitmap(w, h);
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
            DFTImage->SetPixel(j, i, Color::FromArgb(valuePixeli, valuePixeli, valuePixeli));
        }
    }
    pictureBox_OutputImage->Image = DFTImage;
}

System::Void MyForm::inverseDiscreteFourierTransformToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
    int w = dataManager->GetImageWidth();
    int h = dataManager->GetImageHeight();

    // 利用傅立葉之平移性，平移頻率
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            int valuePixeli = dataManager->GetInputImage()[i][j];
            valuePixeli = valuePixeli * pow((float)-1, (float)(i + j));
            dataManager->SetPixel(j, i, valuePixeli);
        }
    }
    fourierTransformMethod->InverseDiscreteFourierTransform(dataManager->GetInputImage(), dataManager->GetOutputImage(), dataManager->GetFreqReal(), dataManager->GetFreqImag(), h, w);

    //將算出頻率資訊傳入輸出影像
    Bitmap^ IDFTImage = gcnew Bitmap(w, h);
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
            IDFTImage->SetPixel(j, i, Color::FromArgb(valuePixeli, valuePixeli, valuePixeli));
        }
    }

    pictureBox_OutputImage->Image = IDFTImage;
}

System::Void MyForm::setResultImageAsSourceImageToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
    int w = dataManager->GetImageWidth();
    int h = dataManager->GetImageHeight();
    Bitmap^ sImage = gcnew Bitmap(w, h);

    //將當前輸出影像作為輸入影像
    for (int i = 0; i < w; i++)
    {
        for (int j = 0; j < h; j++)
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
            dataManager->SetPixel(j, i, valuePixeli);
            sImage->SetPixel(j, i, Color::FromArgb(valuePixeli, valuePixeli, valuePixeli));
        }
    }
    pictureBox_SourceImage->Image = sImage;
}

[STAThread]
int main(array<String^>^ argv)
{
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false);
	FourierTransform::MyForm windowsForm;
	Application::Run(%windowsForm);
}
