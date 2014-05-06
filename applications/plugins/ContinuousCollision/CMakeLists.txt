include(${SOFA_CMAKE_DIR}/preProject.cmake)

set(HEADER_FILES

    PolytopeModel.h
    PolytopeModel.hpp
    initContinuousCollision.h

	)

set(SOURCE_FILES

    initContinuousCollision.cpp
    PolytopeModel.cpp

	)

AddCompilerDefinitions("SOFA_BUILD_PLUGIN_CONTINUOUS_COLLISION")

set_source_files_properties(${SOURCE_FILES} PROPERTIES COMPILE_FLAGS -std=c++11)
add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
include(${SOFA_CMAKE_DIR}/postProject.cmake)
