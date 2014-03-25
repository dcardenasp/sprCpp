/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "vector"

const    unsigned int    Dimension = 3;
typedef  float           IntensitySpelType;
typedef  float           LabelSpelType;
const double PI = std::atan(1.0)*4;

//  The types of the input images are instantiated by the following lines.
//
typedef itk::Image< LabelSpelType, Dimension >              LabelImageType;

typedef itk::ImageFileReader< LabelImageType  >             LabelImageReaderType;
typedef itk::ImageFileWriter< LabelImageType  >             LabelImageWriterType;

void CreateLabelImage(LabelImageType::Pointer p, LabelImageType::Pointer refImg)
{
  p->SetRegions(refImg->GetLargestPossibleRegion());
  p->SetOrigin(refImg->GetOrigin());
  p->SetSpacing(refImg->GetSpacing());
  p->SetDirection(refImg->GetDirection());
  p->Allocate();
  p->FillBuffer(static_cast<LabelImageType::PixelType>(0));
}


int main( int argc, char *argv[] )
{
  if( argc < 4 )
  {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage:\n" << argv[0];
    std::cerr << " numberOfLabels listOfSegmentations output_label_image";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  unsigned int C = std::atoi(argv[1]);

  std::vector<LabelImageType::Pointer> images;


  /*
  //Label template reader
  LabelImageReaderType::Pointer lblReader   = IntensityImageReaderType::New();
  lblReader->SetFileName(argv[2]);

  //Query image reader
  IntensityImageReaderType::Pointer imgReader   = IntensityImageReaderType::New();
  imgReader->SetFileName(argv[3]);

  try
  {
    tptReader->Update();
    lblReader->Update();
    imgReader->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Exception thrown while reading the images" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }

  //Intensity normalization
  MRNormalizationFilterType::Pointer tptNormalizer = MRNormalizationFilterType::New();
  tptNormalizer->SetInput( tptReader->GetOutput() );
  tptNormalizer->Update();

  MRNormalizationFilterType::Pointer imgNormalizer = MRNormalizationFilterType::New();
  imgNormalizer->SetInput( imgReader->GetOutput() );
  imgNormalizer->Update();


  // Create components required for registration
  MetricType::Pointer                   metric          = MetricType::New();
  TransformType::Pointer                transform       = TransformType::New();
  OptimizerType::Pointer                optimizer       = OptimizerType::New();
  IntensityInterpolatorType::Pointer    interpolator    = IntensityInterpolatorType::New();
  RegistrationType::Pointer             registration    = RegistrationType::New();

  //  Initialize the affine transform by using Image Intensity Moments
  TransformInitializerType::Pointer     initializer     = TransformInitializerType::New();
  initializer->SetTransform(   transform );
  initializer->SetFixedImage(  imgNormalizer->GetOutput() );
  initializer->SetMovingImage( tptNormalizer->GetOutput() );
  initializer->MomentsOn();
  initializer->InitializeTransform();
  registration->SetInitialTransformParameters( transform->GetParameters() );

  // Setup the Mattes mutual information metric
  metric->SetNumberOfHistogramBins( 100 );
  const unsigned int numberOfPixels  = imgNormalizer->GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels();
  const unsigned int numberOfSamples = static_cast< unsigned int >( numberOfPixels * 0.1 );
  metric->SetNumberOfSpatialSamples( numberOfSamples );

  // Each component is now connected to the instance of the registration method.
  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetTransform(     transform     );
  registration->SetInterpolator(  interpolator  );

  // Set the registration inputs
  registration->SetFixedImage(       imgNormalizer->GetOutput() );
  registration->SetMovingImage(      tptNormalizer->GetOutput() );
  registration->SetFixedImageRegion( imgNormalizer->GetOutput()->GetLargestPossibleRegion() );

  optimizer->SetMaximumStepLength( .1 ); // If this is set too high, you will get a
  //"itk::ERROR: MeanSquaresImageToImageMetric(0xa27ce70): Too many samples map outside moving image buffer: 1818 / 10000" error

  optimizer->SetMinimumStepLength( 0.01 );

  // Set a stopping criterion
  optimizer->SetNumberOfIterations( 200 );

  try
  {
    registration->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }

  LabelInterpolatorType::Pointer    lblInterpolator = LabelInterpolatorType::New();
  LabelResampleFilterType::Pointer  lblResampler    = LabelResampleFilterType::New();
  lblResampler->SetInput(           lblReader->GetOutput() );
  lblResampler->SetInterpolator(    lblInterpolator );
  lblResampler->SetTransform(       registration->GetOutput()->Get() );
  lblResampler->SetSize(            imgNormalizer->GetOutput()->GetLargestPossibleRegion().GetSize() );
  lblResampler->SetOutputOrigin(    imgNormalizer->GetOutput()->GetOrigin() );
  lblResampler->SetOutputSpacing(   imgNormalizer->GetOutput()->GetSpacing() );
  lblResampler->SetOutputDirection( imgNormalizer->GetOutput()->GetDirection() );  
  lblResampler->SetDefaultPixelValue( 0 );

  LabelImageWriterType::Pointer lblWriter = LabelImageWriterType::New();
  lblWriter->SetInput( lblResampler->GetOutput() );
  lblWriter->SetFileName(argv[4]);

  try
  {
    lblWriter->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Exception thrown while writing the image" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }*/

  return EXIT_SUCCESS;
}
