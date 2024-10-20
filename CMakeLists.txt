# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Set the project name
project(MyCppProject)

# Set C++ standard (optional)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the source and include directories
set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

# Include the directory for header files
include_directories(${INCLUDE_DIR})

# Set the output directory for executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Create the bin directory if it doesn't exist
file(MAKE_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# Use ExternalProject to download and build jemalloc if it isn't found
include(ExternalProject)

# Try to find jemalloc first on the system
find_library(JEMALLOC_LIB jemalloc)

# If jemalloc is not found, download and build it using ExternalProject
if(NOT JEMALLOC_LIB)
    message(STATUS "jemalloc not found, downloading and building it...")

    ExternalProject_Add(
        jemalloc
        GIT_REPOSITORY https://github.com/jemalloc/jemalloc.git
        GIT_TAG 5.3.0  # Specify the version you want to use
        CONFIGURE_COMMAND ./autogen.sh --prefix=${CMAKE_BINARY_DIR}/jemalloc_install --with-private-libpthread
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        BUILD_IN_SOURCE 1
        PREFIX ${CMAKE_BINARY_DIR}/_deps/jemalloc
    )

    # Set the path to jemalloc built library and include directories
    ExternalProject_Get_Property(jemalloc install_dir)
    set(JEMALLOC_LIB "${CMAKE_BINARY_DIR}/jemalloc_install/lib/libjemalloc.a")
    include_directories("${CMAKE_BINARY_DIR}/jemalloc_install/include")

    # Create a target for jemalloc
    add_custom_target(libjemalloc DEPENDS ${JEMALLOC_LIB})

    # Set jemalloc library path
    set(JEMALLOC_LIBS ${JEMALLOC_LIB} pthread dl)  # Added dl for dlsym
else()
    set(JEMALLOC_LIBS ${JEMALLOC_LIB})
endif()

# Find all .cpp files in the src directory
file(GLOB CPP_FILES "${SRC_DIR}/e-centric/*.cpp")

# Loop through each .cpp file and create an executable
foreach(CPP_FILE ${CPP_FILES})
    # Get the filename without the extension
    get_filename_component(EXE_NAME ${CPP_FILE} NAME_WE)
    
    # Add executable for each source file
    add_executable(${EXE_NAME} ${CPP_FILE})
    
    # Ensure the executable is placed in the bin directory
    set_target_properties(${EXE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
    
    # Link jemalloc (with pthread and dl internally) but no pthread for executable
    if(NOT JEMALLOC_LIB)
        add_dependencies(${EXE_NAME} jemalloc)  # Ensure jemalloc is built before executables
    endif()
    
    # Ensure executables depend on jemalloc library if needed
    add_dependencies(${EXE_NAME} libjemalloc)
    
    # Link jemalloc, pthread, and dl to the executable
    target_link_libraries(${EXE_NAME} PRIVATE ${JEMALLOC_LIBS})
endforeach()

