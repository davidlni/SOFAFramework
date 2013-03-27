cmake_minimum_required(VERSION 2.8)

project("SofaGuiCommon")

include(${SOFA_CMAKE_DIR}/pre.cmake)

set(HEADER_FILES

	BaseGUI.h
    BaseViewer.h
	BatchGUI.h
	ColourPickingVisitor.h
	MouseOperations.h
	OperationFactory.h
	PickHandler.h
	FilesRecentlyOpenedManager.h
	SofaGui.h
	ViewerFactory.h
	GUIManager.h
	)

set(SOURCE_FILES

	BaseGUI.cpp
    BaseViewer.cpp
	BatchGUI.cpp
	ColourPickingVisitor.cpp
	FilesRecentlyOpenedManager.cpp
	MouseOperations.cpp
	PickHandler.cpp
	GUIManager.cpp
	ViewerFactory.cpp
	)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})

set(COMPILER_DEFINES "SOFA_BUILD_SOFAGUI")
set(LINKER_DEPENDENCIES SofaComponentMain SofaGraphComponent SofaBaseCollision SofaUserInteraction SofaBaseVisual)

include(${SOFA_CMAKE_DIR}/post.cmake)
