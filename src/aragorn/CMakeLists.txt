find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

add_custom_command(
    OUTPUT 
        aragorn.c
        aragorn.h
        data.c
    COMMENT
        "Converting ${ARAGORN_C} to a library structure"
    COMMAND
        ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/scripts/libraryfier.py 
        -i ${ARAGORN_C}
        -o ${CMAKE_CURRENT_BINARY_DIR}
        --cpp ${CMAKE_C_COMPILER}
    MAIN_DEPENDENCY
        ${ARAGORN_C}
    DEPENDS
        ${PROJECT_SOURCE_DIR}/src/scripts/libraryfier.py
    VERBATIM
)

set(LIBARAGORN_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/aragorn.c
    ${CMAKE_CURRENT_BINARY_DIR}/aragorn.h
    ${CMAKE_CURRENT_BINARY_DIR}/data.c
)

add_library(aragorn ${LIBARAGORN_SOURCES})
target_include_directories(aragorn PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
