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
    //�N�v��Ū�J�ó]�x�s��DataManager����Ƶ��c��
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
            Color srcColor = OriginalImage->GetPixel(j, i); // �^���C���I���C��
            int srcGrey = srcColor.R*0.299 + srcColor.G*0.587 + srcColor.B*0.144; // �m��T�q�D�ন�Ƕ�
            dataManager->SetPixel(j, i, srcGrey);
        }
    }

    std::cout << "-Image has been loaded-" << std::endl;
}

System::Void MyForm::discreteFourierTransformToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
    int w = dataManager->GetImageWidth();
    int h = dataManager->GetImageHeight();

    // �Q�γť߸��������ʡA�����W�v
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            int valuePixeli = dataManager->GetInputImage()[i][j];
            valuePixeli = valuePixeli * pow((float)-1, (float)(i + j));
            dataManager->SetPixel(j, i, valuePixeli);
        }
    }

    //�N��X�W�v��T�ǤJ��X�v��
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

    // �Q�γť߸��������ʡA�����W�v
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

    //�N��X�W�v��T�ǤJ��X�v��
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

    //�N��e��X�v���@����J�v��
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
