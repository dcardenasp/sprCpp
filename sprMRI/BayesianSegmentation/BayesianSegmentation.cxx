#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkComposeImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include "itkCenteredTransformInitializer.h"
#include "itkVersorRigid3DTransform.h"
#include "itkAffineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkImageRegistrationMethod.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"

#include "itkImageDuplicator.h"
#include "itkVectorImageToImageAdaptor.h"
#include "itkStatisticsImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkBayesianClassifierImageFilter.h"

#include <math.h>
#define PI 3.14159265

//#include "itkBayesianClassifierImageFilter.h"


int main(int argc, char *argv[])
{
  if (argc < 2)
    {
    std::cerr << "Usage: " << argv[0] << " imageFilename priorList" << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImageFileName         = argv[1];
  const char * inputMovingImageFileNames  = argv[2];
  const char * outputImageFileName        = argv[3];
  
  typedef itk::Image<float, 3>  ScalarImageType;
  typedef itk::Image<unsigned char, 3>  LabelImageType;

  typedef itk::ImageFileReader<ScalarImageType> ReaderType;
  ReaderType::Pointer fixed = ReaderType::New();
  fixed->SetFileName(inputImageFileName);
  fixed->Update();

  std::cout << "Input read\n";

  //*******LOAD AND CREATE PRIOR VECTOR IMAGE *********//
  std::vector<ScalarImageType::Pointer> prImages;
  std::ifstream infile(inputMovingImageFileNames);
  std::string line;
  int numClasses = 0;
  std::getline(infile, line);
  ReaderType::Pointer moving = ReaderType::New();
  moving->SetFileName(line);
  moving->Update();

  while (std::getline(infile, line))
  {
      ReaderType::Pointer tmp = ReaderType::New();
      tmp->SetFileName(line);
      tmp->Update();

      prImages.push_back(tmp->GetOutput());
      numClasses++;
  }
  std::cout << "Prior read and built\n";

  //******AFFINE REGISTERING PRIORS**********//
  typedef itk::VersorRigid3DTransform< double > RigidTransformType;
  typedef itk::AffineTransform< double, 3 > AffineTransformType;
  typedef itk::CenteredTransformInitializer< RigidTransformType,
          ScalarImageType, ScalarImageType > TransformInitializerType;
  typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
  typedef itk::MattesMutualInformationImageToImageMetric<
          ScalarImageType, ScalarImageType > MetricType;

  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  RigidTransformType::Pointer  rigidTransform = RigidTransformType::New();
  initializer->SetTransform(   rigidTransform );
  initializer->SetFixedImage(  fixed->GetOutput() );
  initializer->SetMovingImage( moving->GetOutput() );
  initializer->MomentsOn();
  initializer->InitializeTransform();
  std::cout << "Rigid Transform Initialization completed" << std::endl;
  typedef itk:: LinearInterpolateImageFunction<
          ScalarImageType, double > InterpolatorType;
  typedef itk::ImageRegistrationMethod<
          ScalarImageType, ScalarImageType >    RegistrationType;
  MetricType::Pointer               metric        = MetricType::New();
  OptimizerType::Pointer            optimizer     = OptimizerType::New();
  InterpolatorType::Pointer         interpolator  = InterpolatorType::New();
  RegistrationType::Pointer         registration  = RegistrationType::New();
  registration->SetMetric(          metric        );
  registration->SetOptimizer(       optimizer     );
  registration->SetInterpolator(    interpolator  );
  registration->SetFixedImage(      fixed->GetOutput()  );
  registration->SetMovingImage(     moving->GetOutput() );
  metric->SetNumberOfHistogramBins( 50 );
  ScalarImageType::RegionType fixedRegion = fixed->GetOutput()->GetBufferedRegion();
  const unsigned int numberOfPixels = fixedRegion.GetNumberOfPixels();
  metric->ReinitializeSeed( 76926294 );
  registration->SetFixedImageRegion( fixedRegion );
  registration->SetInitialTransformParameters( rigidTransform->GetParameters() );
  registration->SetTransform( rigidTransform );
  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( rigidTransform->GetNumberOfParameters() );
  const double translationScale = 1.0 / 1000.0;
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = translationScale;
  optimizerScales[4] = translationScale;
  optimizerScales[5] = translationScale;
  optimizer->SetScales( optimizerScales );
  optimizer->SetMaximumStepLength( 0.2000  );
  optimizer->SetMinimumStepLength( 0.0001 );
  optimizer->SetNumberOfIterations( 200 );
  metric->SetNumberOfSpatialSamples( 10000L );
  registration->Update();
  rigidTransform->SetParameters( registration->GetLastTransformParameters() );
  std::cout << "Rigid Transform Completed: "
          << registration->GetOptimizer()->GetStopConditionDescription()
          << std::endl;

  AffineTransformType::Pointer  affineTransform = AffineTransformType::New();
  affineTransform->SetCenter( rigidTransform->GetCenter() );
  affineTransform->SetTranslation( rigidTransform->GetTranslation() );
  affineTransform->SetMatrix( rigidTransform->GetMatrix() );
  registration->SetTransform( affineTransform );
  registration->SetInitialTransformParameters( affineTransform->GetParameters() );
  optimizerScales = OptimizerScalesType( affineTransform->GetNumberOfParameters() );
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = 1.0;
  optimizerScales[4] = 1.0;
  optimizerScales[5] = 1.0;
  optimizerScales[6] = 1.0;
  optimizerScales[7] = 1.0;
  optimizerScales[8] = 1.0;
  optimizerScales[9]  = translationScale;
  optimizerScales[10] = translationScale;
  optimizerScales[11] = translationScale;
  optimizer->SetScales( optimizerScales );
  optimizer->SetMaximumStepLength( 0.2000  );
  optimizer->SetMinimumStepLength( 0.0001 );
  optimizer->SetNumberOfIterations( 200 );
  metric->SetNumberOfSpatialSamples( 50000L );
  registration->Update();
  affineTransform->SetParameters( registration->GetLastTransformParameters() );
  std::cout << "Affine Registration completed" << std::endl;
  //typedef itk::Image<itk::Vector< float, 6>, 3> VectorImageType;
  typedef itk::ComposeImageFilter<ScalarImageType> ImageToVectorImageFilterType;
  ImageToVectorImageFilterType::Pointer vectorBuilder = ImageToVectorImageFilterType::New();
  //std::vector<ScalarImageType::Pointer> posImages;
  typedef itk::ResampleImageFilter< ScalarImageType,
          ScalarImageType>    ResampleFilterType;
  ResampleFilterType::Pointer resample = ResampleFilterType::New();
  resample->SetTransform( affineTransform );
  resample->SetSize(    fixed->GetOutput()->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin(  fixed->GetOutput()->GetOrigin() );
  resample->SetOutputSpacing( fixed->GetOutput()->GetSpacing() );
  resample->SetOutputDirection( fixed->GetOutput()->GetDirection() );
  resample->SetInterpolator( interpolator );
  resample->SetDefaultPixelValue(0.0);
  for(unsigned int c=0; c<numClasses; c++)
  {
      resample->SetInput( prImages[c] );
      resample->Update();
      prImages[c] = resample->GetOutput();
      prImages[c]->DisconnectPipeline();
      vectorBuilder->SetInput(c,prImages[c]);
      resample->SetDefaultPixelValue(0.0);
  }
  vectorBuilder->Update();
  typedef ImageToVectorImageFilterType::OutputImageType VectorImageType;
  VectorImageType::Pointer priors = vectorBuilder->GetOutput();
  typedef itk::ImageDuplicator< VectorImageType > DuplicatorType;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(vectorBuilder->GetOutput());
  duplicator->Update();
  VectorImageType::Pointer posteriors = duplicator->GetOutput();

  //*******Partitioning the image********//
  itk::ImageRegionConstIterator<ScalarImageType> itFixed(fixed->GetOutput(),fixed->GetOutput()->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<VectorImageType> itPr(priors,priors->GetLargestPossibleRegion());
  itk::ImageRegionIterator<VectorImageType> itPos(posteriors,posteriors->GetLargestPossibleRegion());
  VectorImageType::PixelType sumPr, h, means, cov, tmp;
  sumPr.SetSize(numClasses);
  h.SetSize(numClasses);
  means.SetSize(numClasses);
  cov.SetSize(numClasses);
  tmp.SetSize(numClasses);
  sumPr.Fill(0.0);

  for(itPr.GoToBegin(); !itPr.IsAtEnd(); ++itPr)
      sumPr = sumPr + itPr.Get();

  for(unsigned int iter = 0; iter<atoi(argv[4]); iter++)
  {
    h.Fill(0.0);
    means.Fill(0.0);
    cov.Fill(0.0);
    tmp.Fill(0.0);
    for(itFixed.GoToBegin(), itPos.GoToBegin(); !itFixed.IsAtEnd(); ++itFixed, ++itPos)
    {
        h = h + itPos.Get();
        means = means + itFixed.Get()*itPos.Get();
    }
    float sh=0;
    for(unsigned int c=0; c<numClasses; c++)
    {
        means[c] = means[c]/h[c];
        sh += h[c];
    }
    std::cout << iter << "\t" << h;
    for(itFixed.GoToBegin(), itPos.GoToBegin(); !itFixed.IsAtEnd(); ++itFixed, ++itPos)
    {
        for(unsigned int c=0; c<numClasses; c++)
            tmp[c] = (itPos.Get()[c])*pow(itFixed.Get()-means[c],2.0);
        cov = cov + tmp;
    }
    for(unsigned int c=0; c<numClasses; c++)
        cov[c] = cov[c]/h[c];

    std::cout << "\t" << means << "\t" << cov << std::endl;

    for(itFixed.GoToBegin(), itPr.GoToBegin(), itPos.GoToBegin(); !itFixed.IsAtEnd(); ++itFixed, ++itPr, ++itPos)
    {
        float sumQ = 0;
        for(unsigned int c=0; c<numClasses; c++)
        {
            tmp[c] = itFixed.Get();
            tmp[c] = pow(2*PI*cov[c],-0.5)*exp(-pow(tmp[c]-means[c],2.0)/(2*cov[c]));
            //tmp[c] = tmp[c]*h[c];
            if(itPr.Get()[c]>0)
                tmp[c] = tmp[c]*h[c]*(itPr.Get()[c])/sumPr[c];
            else
                tmp[c] = tmp[c]*h[c]*(1e-4)/sumPr[c];
            sumQ += tmp[c];
        }

        for(unsigned int c=0; c<numClasses; c++)
            tmp[c] = tmp[c]/sumQ;
        itPos.Set(tmp);
    }
  }
  std::cout << "Partitioning completed" << std::endl;

  //*******Labeling the image********//
  LabelImageType::Pointer lbl = LabelImageType::New();
  lbl->SetRegions(   fixed->GetOutput()->GetLargestPossibleRegion() );
  lbl->SetOrigin(    fixed->GetOutput()->GetOrigin() );
  lbl->SetSpacing(   fixed->GetOutput()->GetSpacing() );
  lbl->SetDirection( fixed->GetOutput()->GetDirection() );
  lbl->Allocate();
  itk::ImageRegionIterator<LabelImageType> itLbl(lbl,lbl->GetLargestPossibleRegion());
  for(itLbl.GoToBegin(), itPos.GoToBegin(); !itPos.IsAtEnd(); ++itLbl, ++itPos)
  {
      unsigned char etiq = 0;
      for(unsigned int c=1; c<numClasses; c++)
      {
          etiq = (itPos.Get()[c]>itPos.Get()[etiq])?c:etiq;
      }
      itLbl.Set(etiq);
  }
  //***************WRITER**************//

  typedef itk::ImageFileWriter< LabelImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(lbl);
  writer->SetFileName(outputImageFileName);
  writer->Update();

  //std::cout << pr->GetNumberOfComponentsPerPixel() << std::endl;
  //std::cout << reader->GetOutput()->GetLargestPossibleRegion() << std::endl;
  
  return EXIT_SUCCESS;
}
