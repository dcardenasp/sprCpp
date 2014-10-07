CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(BayesianSegmentation)
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

include_directories(../MRIntensityNormalization)

add_executable(bayes_train MACOSX_BUNDLE sprBayesTrain.cxx )
add_executable(BayesianSegmentation MACOSX_BUNDLE BayesianSegmentation.cxx)

if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(BayesianSegmentation ITKReview ${ITK_LIBRARIES})
  target_link_libraries(bayes_train ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(BayesianSegmentation ${ITK_LIBRARIES})
  target_link_libraries(bayes_train ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )
