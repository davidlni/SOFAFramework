include(${SOFA_CMAKE_DIR}/preProject.cmake)

set(HEADER_FILES

    initContinuousCollision.h

	)

set(SOURCE_FILES

    initContinuousCollision.cpp

	)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})

include(${SOFA_CMAKE_DIR}/postProject.cmake)
