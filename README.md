# GraCFL

## Table of Contents
- [Project Overview](#project-overview)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Build Instructions](#build-instructions)
- [Run Scripts](#run-scripts)
- [Benchmark Graphs and Grammar Files](#benchmark-graphs-and-grammar-files)
- [Example Structure](#example-structure)
- [Run with Singularity](#run-with-singularity)

## Project Overview
**GraCFL** is designed for high-performance CFL reachability analysis. The project includes various executables in the `build/bin/` directory that perform CFL reachability computations with different optimizations and parallel strategies.

## Dependencies
The project depends on the following libraries and tools:

- **CMake** (minimum version 3.26)
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
    For parallel executables, run the following command
   ```bash
   ./<executable_name> <graph_file> <grammar_file> <no_threads>
   ```
   
    - **`<executable_name>`**: The name of the executable you wish to run (e.g., `e-centric-bi-grammar_driven-parallel`).
    - **`<graph_file>`**: Path to the input graph file.
    - **`<grammar_file>`**: Path to the input grammar file.
    - **`<no_threads>`**: No of threads to run in parallel (only applicable for the parallel executables)


## Run Scripts 

The project includes scripts to help automate the process of running the generated executables on specific input files. These scripts are designed to make it easier to reproduce the results used in the related research paper or run all available executables. 

The scripts are located inside the `GraCFL/` directory. Run these scripts only after building the project.

### Running the Executables Used in the Paper

To generate results for the executables that were specifically used in the research paper, you can run the `run_script_paper.sh` file. This script will automatically run the relevant executables with their corresponding input files (graphs and grammars - will have to be put in place manually).

To execute the scripts:

```bash
sh run_script_paper.sh
```
The generated results may differ from those reported in the paper due to variations in your system configuration.

### Running All the Executables

You can run the `run_script_all.sh` file to generate results for all the executables. `run_script_all.sh` runs all executables present in the `/build/bin/` directory, even those that were not used in the paper.

To execute the script:

```bash
sh run_script_all.sh
```

## Benchmark Graphs and Grammar Files

The benchmark graphs and grammar files used in the research paper are available for download. These files are necessary to reproduce the experiments and results described in the paper. You can download them from the following Google Drive link:

[Download Benchmark Graphs and Grammar Files](https://drive.google.com/drive/folders/1mLMhuXbj4Eu9xrhFcQs3uH3Lp0lw0iZ8?usp=sharing)

### How to Use

1. Download the graph and grammar files from the provided Google Drive link.
2. Place the downloaded `graphs_grammars` folder in the `GraCFL` directory:
   - Graph files are located inside  the `graphs_grammars/graphs/` directory.
   - Grammar files are located inside the `graphs_grammars/grammars/` directory.

Once the files are in place, you can run the provided scripts (e.g., `run_script_paper.sh` or `run_script_all.sh`) to generate results based on these input files.

Make sure the directory structure matches the paths used in the scripts to ensure proper execution.

## Example Structure

Below is the example directory structure of the project:
```
GraCFL/
├── include/                # Header files
├── src/                    # Source files
├── graphs_grammars/        # Graph and Grammar input files
├── build/bin/              # Generated executables (after build)
├── logs/                   # Log files (generated during runtime)
├── CMakeLists.txt          # CMake build configuration
├── README.md               # This file
├── run_script_paper.sh     # Script to automate the execution of the executables used in the paper
└── run_script_all.sh       # Script to automate the execution of all executables
```

## Run with Singularity

You can run `GraCFL` using the provided Singularity image without needing to build it on your machine.

[Download Singularity Image](https://drive.google.com/file/d/1n9EJr4GsLLwcBRBuU-gkIfkpyyDJwc3N/view?usp=sharing) 

### How to Use

1. Download the Singularity image
2. Make sure you have `singularity` installed on your machine
3. Download the `graphs_grammars` directory from the provided link.
4. Run the following command to start the container:
```bash
   singularity shell \
   --bind <path_to_graphs_grammars_parent_dir>/graphs_grammars/:/usr/src/mycppproject/graphs_grammars/ \
   --bind <path_to_graphs_grammars_parent_dir>/logs/:/usr/src/mycppproject/logs/ \
   gracfl_amd64_1_0.sif
```
Here:
- `<path_to_graphs_grammars_parent_dir>/graphs_grammars/` should be the path to your `graphs_grammars` directory.
- `<path_to_graphs_grammars_parent_dir>/logs/` should be the path to a directory where you want to store the output. Make sure to create this `logs` directory in advance by running:
```bash 
   mkdir -p logs
```
5. Now run the following commands to run the `run_script_paper.sh` or `run_script_all.sh`
```bash
   cd /usr/src/mycppproject/   # Navigate to the working directory
   ./run_script_paper.sh       # Use this to reproduce the results from the paper
   ./run_script_all.sh         # Use this to run all methods (including those not covered in the paper)
```
6. The output logs will be available inside the `logs` directory you specified.






