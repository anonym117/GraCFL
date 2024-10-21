# GraCFL

## Table of Contents
- [Project Overview](#project-overview)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Build Instructions](#build-instructions)
- [Usage](#usage)
- [Example](#example)
- [Example Structure](#example-structure)
- [Contributing](#contributing)
- [License](#license)

## Project Overview
**GraCFL** is designed for high-performance CFL reachability analysis. The project includes various executables in the `build/bin/` directory that perform CFL reachability computations with different optimizations and parallel strategies.

## Dependencies
The project depends on the following libraries and tools:

- **CMake** (minimum version 3.10)
- **C++11 or higher** (supported by the project)
- **TBB** (Threading Building Blocks)
- **jemalloc** (for memory allocation)
- **OpenMP** (optional for parallelism)

The project automatically fetches and builds the necessary dependencies if they are not found on the system.

## Installation

### Clone the Repository
To get started, clone the repository using the following command:

```bash
git clone https://github.com/anonym117/GraCFL.git
cd GraCFL
```

## Build Instructions

To build the project, follow these steps:

1. **Create a build directory:**

    First, create a directory for building the project and move into that directory:

    ```bash
    mkdir build
    cd build
    ```

2. **Configure the project using CMake:**

    Run the `cmake` command from the `build` directory to configure the project. By default, the project is set to use Release mode with `-O3` optimization:

    ```bash
    cmake ..
    ```
3. **Build the project:**

    Run the `make` command to build the project:
   
    ```bash
    make
    ```
    This will compile all source files and place the resulting executables in the `bin/` directory under the `build/` directory.

5. **Running the executables:**

    After building, you can run the generated executables from the `build/bin/` directory as follows:
   
   ```bash
   ./<executable_name> <graph_file> <grammar_file>
   ```



