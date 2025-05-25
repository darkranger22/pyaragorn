cmake_minimum_required(VERSION 3.20)
project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES C)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(ARAGORN_VERSION "1.2.41")

# --- Setup Cython -------------------------------------------------------------

set_property(GLOBAL PROPERTY PYTHON_EXTENSIONS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(CYTHON_HEADERS_DIR ${CMAKE_CURRENT_LIST_DIR}/include)
include(src/scripts/cmake/CythonExtension.cmake)

# --- Set ----------------------------------------------------------------------

set(ARAGORN_C ${CMAKE_CURRENT_LIST_DIR}/vendor/aragorn/aragorn${ARAGORN_VERSION}.c)

# --- Compile code -------------------------------------------------------------

add_subdirectory(src)
