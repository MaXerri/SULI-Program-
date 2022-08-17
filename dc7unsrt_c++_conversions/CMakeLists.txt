cmake_minimum_required(VERSION 3.0)

set(DEPENDENCIES
ncdist/B4.cpp          ncdist/B4.h
ncdist/C3.cpp          ncdist/C3.h
ncdist/CellInputData.h ncdist/CellInputData.cpp
ncdist/D7.cpp          ncdist/D7.h
ncdist/D7.h            ncdist/D7.cpp 
ncdist/D7_Boundary.h   ncdist/D7_Boundary.cpp
ncdist/D7_BoundaryList.h ncdist/D7_BoundaryList.cpp 
ncdist/D7_Subboundary.h ncdist/D7_Subboundary.cpp 
ncdist/D7_ClassifySubboundaries.h ncdist/D7_ClassifySubboundaries.cpp
ncdist/G6.cpp  ncdist/G6.h
ncdist/inverse.cpp     ncdist/inverse.h
ncdist/LatticeConverter.h ncdist/LatticeConverter.cpp
ncdist/LRL_Cell.h      ncdist/LRL_Cell.cpp  
ncdist/LRL_Cell_Degrees.h   
ncdist/LRL_Cell_Degrees.cpp  
ncdist/LRL_CoordinateConversionMatrices.h ncdist/LRL_CoordinateConversionMatrices.cpp  
ncdist/LRL_StringTools.h  ncdist/LRL_StringTools.cpp  
ncdist/MatD7.h         ncdist/MatD7.cpp  
ncdist/MatG6.h         ncdist/MatG6.cpp  
ncdist/MatS6.h         ncdist/MatS6.cpp  
ncdist/MatN.cpp        ncdist/MatN.h
ncdist/MatMN.cpp       ncdist/MatMN.h
ncdist/Mat66.cpp       ncdist/Mat66.h
ncdist/MaximaTools.h   ncdist/MaximaTools.cpp
ncdist/PrintTable.h    ncdist/PrintTable.cpp
ncdist/ProjectorTools.h ncdist/ProjectorTools.cpp 
ncdist/RandTools.h     ncdist/RandTools.cpp 
ncdist/ReadCellData.h  ncdist/ReadCellData.cpp
ncdist/Reducer.cpp     ncdist/Reducer.h
ncdist/S6.cpp          ncdist/S6.h
ncdist/S6Dist.cpp      ncdist/S6Dist.h
ncdist/Selling.cpp     ncdist/Selling.h
ncdist/S6M_SellingReduce.h
ncdist/S6Dist_func.cpp ncdist/S6Dist_func.h
ncdist/Vec_N_Tools.cpp ncdist/Vec_N_Tools.h
ncdist/VecN.cpp        ncdist/VecN.h
ncdist/VectorTools.cpp ncdist/VectorTools.h
ncdist/vector_3d.cpp   ncdist/vector_3d.h)

set(NCDIST_SOURCES
ncdist/Delone.h    ncdist/Delone.cpp
${DEPENDENCIES}
)

project(main_executable)

add_executable(main_executable converter_main.cpp conversions.cpp ncdistmatmod.cpp ${NCDIST_SOURCES})
target_include_directories(main_executable PUBLIC ${CMAKE_CURRENT_SOURCE_DIRECTORY}/ncdist PUBLIC ${CMAKE_CURRENT_SOURCE_DIRECTORY})
#target_compile_features(main_executable PRIVATE cxx_std_11)
