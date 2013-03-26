cmake_minimum_required(VERSION 2.8)

project("SofaMiscEngine")

include(${SOFA_CMAKE_DIR}/pre.cmake)

set(HEADER_FILES

    initMiscEngine.h 
    engine/Distances.h 
    engine/Distances.inl

    )
    
set(SOURCE_FILES

    initMiscEngine.cpp 
    engine/Distances.cpp
 
    )
    

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})

set(COMPILER_DEFINES "SOFA_BUILD_MISC_ENGINE" )
set(LINKER_DEPENDENCIES SofaNonUniformFem )
    
include(${SOFA_CMAKE_DIR}/post.cmake)
