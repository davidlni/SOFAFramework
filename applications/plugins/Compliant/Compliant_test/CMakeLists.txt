include(${SOFA_CMAKE_DIR}/preProject.cmake)

set(HEADER_FILES

    Compliant_test.h
	Python_test.h

        )

set(SOURCE_FILES

    Assembly_test.cpp
    AssembledSolver_test.cpp
	Python_test.cpp
	python_test_list.cpp

	)

if(APPLE)
#    set(RC_FILES "runSOFA.icns")
endif()

if(WIN32)
	link_directories("${SOFA_EXTLIBS_DIR}/SuiteSparse/cholmod/Lib")
endif()

add_executable(${PROJECT_NAME} ${HEADER_FILES} ${SOURCE_FILES} ${RC_FILES})

AddLinkerDependencies(Compliant SofaTest)

AddCompilerDefinitions("COMPLIANT_TEST_PYTHON_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/python\"")

# needed by python tests
AddCompilerDefinitions("SOFA_SRC_DIR=${SOFA_SRC_DIR}")
AddLinkerDependencies(SofaPython)


include(${SOFA_CMAKE_DIR}/postProject.cmake)
