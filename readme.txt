!. The singularity image mentioned in the README.md is used for a different g++ --version. The version used in the singularity image is g++ (Ubuntu 5.4.0-6ubuntu1~16.04.9) 5.4.0 20160609
   The command to invoke the singularity image is given below:
        ```
        singularity shell --nv /home/csgrads/sfuad001/graspan-research/singularity-image/ucr-singularity-ml-software-master-latest.simg
        ```
@. The following command compiles the cpp files using the g++ compiler inside singularity image
    ```
    g++ --std=c++11 -fcilkplus -o data-driven data-driven.cpp globals.hpp grammar.hpp 
    ```
    ```
	g++ --std=c++11 -fcilkplus -o topo-driven topo-driven.cpp globals.hpp grammar.hpp
    ```
    ```
	 g++ --std=c++11 -fcilkplus -o ooc  ooc.cpp globals.hpp grammar.hpp
    ```

### log file name convention
type-filename-date-machinename-sequenceno.log
Example: dataflow-httpd-oct12-venus-1.log

# After running the ./data-driven executable, the venus machine got frozen.

This is what I have found so far.


Amir's "data-driven" execution got killed. I think because of the memory shortage. Here is the log:

# Edges in Part 1021: 20202
# Nodes in Part 1022: 10991
# Edges in Part 1022: 19620
# Nodes in Part 1023: 10641
# Edges in Part 1023: 17576
#edges: 18968430 , SumOfParts: 20160606
Killed


