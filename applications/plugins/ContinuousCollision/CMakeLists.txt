include(${SOFA_CMAKE_DIR}/preProject.cmake)

set(HEADER_FILES
    PolytopeModel.h
    PolytopeModel.hpp
    RTriangleModel.h
    RTriangleModel.hpp
    initContinuousCollision.h
    ContinuousDetection.h

	)

set(SOURCE_FILES
    PolytopeModel.cpp
    RTriangleModel.cpp
    initContinuousCollision.cpp
#     ContinuousDetection.cpp
	)

AddCompilerDefinitions("SOFA_BUILD_PLUGIN_CONTINUOUS_COLLISION")

set_source_files_properties(${SOURCE_FILES} PROPERTIES COMPILE_FLAGS -std=c++11)
add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
include(${SOFA_CMAKE_DIR}/postProject.cmake)
