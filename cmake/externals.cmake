cmake_minimum_required(VERSION 2.8)

cmake_policy(SET CMP0015 OLD)

# include dirs
if(WIN32)
    list(APPEND GLOBAL_INCLUDE_DIRECTORIES "${SOFA_INC_DIR}")
endif()
list(APPEND GLOBAL_INCLUDE_DIRECTORIES "${SOFA_FRAMEWORK_DIR}")
list(APPEND GLOBAL_INCLUDE_DIRECTORIES "${SOFA_MODULES_DIR}")
list(APPEND GLOBAL_INCLUDE_DIRECTORIES "${SOFA_APPLICATIONS_DIR}")

if(SOFA-MISC_DEVELOPER_MODE)
	list(APPEND GLOBAL_INCLUDE_DIRECTORIES "${SOFA_APPLICATIONS_DEV_DIR}")
endif()

if(SOFA-EXTERNAL_HAVE_BOOST)
	if(NOT SOFA-EXTERNAL_BOOST_PATH STREQUAL "")
		set(ENV{Boost_DIR} ${SOFA-EXTERNAL_BOOST_PATH})
	endif()
	find_package("Boost" REQUIRED COMPONENTS thread graph system)
	set(Boost_DIR ${SOFA-EXTERNAL_BOOST_PATH} CACHE INTERNAL "Boost root directory" FORCE)
	set(Boost_LIB_DIAGNOSTIC_DEFINITIONS ${Boost_LIB_DIAGNOSTIC_DEFINITIONS} CACHE INTERNAL "Boost lib diagnostic definitions" FORCE)
	list(APPEND GLOBAL_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}")
else()
    list(APPEND GLOBAL_INCLUDE_DIRECTORIES "${SOFA_EXTLIBS_DIR}/miniBoost")
endif(SOFA-EXTERNAL_HAVE_BOOST)

if(SOFA-EXTERNAL_HAVE_EIGEN2)
	list(APPEND GLOBAL_INCLUDE_DIRECTORIES "${SOFA_EXTLIBS_DIR}/eigen-3.1.1")
endif()
list(APPEND GLOBAL_INCLUDE_DIRECTORIES "${SOFA_EXTLIBS_DIR}/newmat")

if(SOFA-EXTERNAL_HAVE_FLOWVR)
	list(APPEND GLOBAL_INCLUDE_DIRECTORIES "${SOFA_EXTLIBS_DIR}/miniFlowVR/include")
endif()

## Zlib (SOFA-EXTERNAL_HAVE_ZLIB)
if(WIN32 OR XBOX)
	set(ZLIB_LIBRARIES "zlib")
else()
	find_library(ZLIB_LIBRARIES "z")
endif()
set(ZLIB_LIBRARIES ${ZLIB_LIBRARIES} CACHE INTERNAL "ZLib Library")
if (SOFA-EXTERNAL_HAVE_ZLIB)
  set(ZLIB_LIBRARIES_OPTIONAL ${ZLIB_LIBRARIES} CACHE INTERNAL "ZLib Library")
else()
  set(ZLIB_LIBRARIES_OPTIONAL "" CACHE INTERNAL "ZLib Library (optional)")
endif()
RegisterDependencies(${ZLIB_LIBRARIES} OPTION SOFA-EXTERNAL_HAVE_ZLIB COMPILE_DEFINITIONS SOFA_HAVE_ZLIB)
# lib dir
link_directories("${SOFA_LIB_DIR}")
link_directories("${SOFA_LIB_OS_DIR}")

# packages and libraries

## opengl / glew / glut
if (NOT SOFA-MISC_NO_OPENGL)
	find_package(OPENGL REQUIRED)
	if(WIN32)
		set(OPENGL_LIBRARIES ${OPENGL_LIBRARIES} "glu32")
		set(GLEW_LIBRARIES "glew32")
		set(GLUT_LIBRARIES "glut32")
		set(PNG_LIBRARIES "libpng")
	else()
		find_package(GLEW REQUIRED)
		find_library(GLUT_LIBRARIES "glut")
		if(SOFA-EXTERNAL_PNG_SPECIFIC_VERSION)
			set(PNG_LIBRARIES "${SOFA-EXTERNAL_PNG_VERSION}")
		else()
			find_library(PNG_LIBRARIES "png")
		endif()
	endif()

  ## GLU
	if(UNIX)
		if(NOT APPLE)
			list(APPEND GLUT_LIBRARIES GLU X11)
		endif()
    list(APPEND GLUT_LIBRARIES dl)
	endif()
	list(REMOVE_DUPLICATES GLUT_LIBRARIES)
else()
	set(OPENGL_LIBRARIES "")
	set(GLEW_LIBRARIES "")
	set(GLUT_LIBRARIES "")
	set(PNG_LIBRARIES "")
endif()

set(OPENGL_LIBRARIES ${OPENGL_LIBRARIES} CACHE INTERNAL "OpenGL Library")
set(GLEW_LIBRARIES ${GLEW_LIBRARIES} CACHE INTERNAL "GLEW Library")
set(GLUT_LIBRARIES ${GLUT_LIBRARIES} CACHE INTERNAL "GLUT Library")
set(PNG_LIBRARIES ${PNG_LIBRARIES} CACHE INTERNAL "PNG Library")

RegisterDependencies(${OPENGL_LIBRARIES})
RegisterDependencies(${GLEW_LIBRARIES} OPTION SOFA-EXTERNAL_HAVE_GLEW COMPILE_DEFINITIONS SOFA_HAVE_GLEW)
RegisterDependencies(${GLUT_LIBRARIES})
RegisterDependencies(${PNG_LIBRARIES} OPTION SOFA-EXTERNAL_HAVE_PNG COMPILE_DEFINITIONS SOFA_HAVE_PNG)

# unit tests
if(SOFA-MISC_TESTS)
	enable_testing()
endif()
