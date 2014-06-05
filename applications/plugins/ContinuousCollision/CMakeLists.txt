include(${SOFA_CMAKE_DIR}/preProject.cmake)

set(HEADER_FILES
    PolytopeModel.h
    PolytopeModel.hpp
    RTriangleModel.h
    RTriangleModel.hpp
    initContinuousCollision.h
    ContinuousDetection.h
    ContinuousIntersection.h
    ContinuousContactMapper.h
    ContinuousPenalityContact.h
	)

set(SOURCE_FILES
    PolytopeModel.cpp
    RTriangleModel.cpp
    initContinuousCollision.cpp
    ContinuousDetection.cpp
    ContinuousIntersection.cpp
    ContinuousContactMapper.cpp
    ContinuousPenalityContact.cpp
	)

AddCompilerDefinitions("SOFA_BUILD_PLUGIN_CONTINUOUS_COLLISION")
AddLinkerDependencies( SofaComponentMisc SofaHelper SofaSimulationGraph SofaGuiMain)

set_source_files_properties(${SOURCE_FILES} PROPERTIES COMPILE_FLAGS -std=c++11)
add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
include(${SOFA_CMAKE_DIR}/postProject.cmake)
