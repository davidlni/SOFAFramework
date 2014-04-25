include(${SOFA_CMAKE_DIR}/preProject.cmake)

set(SOURCE_FILES
#   PolytopeModel_test.cpp
#     ContinuousCollision_test.cpp
  test_CollisionPipeline.cpp
	)

if(APPLE)
    set(RC_FILES "runSOFA.icns")
endif()

# include_directories("${SOFA_EXTLIBS_DIR}/gtest/include")
# link_directories("${SOFA_EXTLIBS_DIR}/gtest/lib")

if(WIN32)
	link_directories("${SOFA_EXTLIBS_DIR}/SuiteSparse/cholmod/Lib")
endif()

include_directories(/home/rortiz/projects/Sofa/extlibs/tinyxml)
set_source_files_properties(${SOURCE_FILES} PROPERTIES COMPILE_FLAGS -std=c++11)
add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${RC_FILES})

# AddLinkerDependencies(ContinuousCollision SofaComponentMain SofaSimulationGraph SofaTest)
AddLinkerDependencies(ContinuousCollision SofaComponentMain SofaSimulationGraph tinyxml)

AddCompilerDefinitions("CCD_TEST_SCENES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/scenes\"")


include(${SOFA_CMAKE_DIR}/postProject.cmake)
