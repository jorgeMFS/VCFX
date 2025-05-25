cmake_minimum_required(VERSION 3.14)

if(NOT PYTHON_BINDINGS)
    return()
endif()

find_package(Python3 COMPONENTS Development REQUIRED)

add_library(_vcfx MODULE bindings.cpp)
target_link_libraries(_vcfx PRIVATE vcfx_core Python3::Python)

# Place the compiled module into the build/python directory
set_target_properties(_vcfx PROPERTIES
    PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python/vcfx"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python/vcfx"
)

configure_file(__init__.py "${CMAKE_BINARY_DIR}/python/vcfx/__init__.py" COPYONLY)
configure_file(tools.py "${CMAKE_BINARY_DIR}/python/vcfx/tools.py" COPYONLY)

install(TARGETS _vcfx
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/vcfx
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/vcfx)
install(FILES __init__.py tools.py DESTINATION ${CMAKE_INSTALL_LIBDIR}/vcfx)
