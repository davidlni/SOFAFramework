include(${SOFA_CMAKE_DIR}/preProject.cmake)

set(SOURCE_FILES

    mixedPendulum.cpp

    )

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

AddLinkerDependencies(SofaSimulationTree SofaComponentMain SofaGuiMain)

include(${SOFA_CMAKE_DIR}/postProject.cmake)

