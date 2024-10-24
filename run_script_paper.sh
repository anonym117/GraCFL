#!/bin/bash

### This script runs a subset of executables to generate results as presented in the paper.

# Set the paths for the executables and graphs directories, and grammar files
BIN_DIR="./build/bin"
GRAPH_DIR="./graphs_grammars/graphs"
GRAPHS_DIR_DATAFLOW="$GRAPH_DIR/dataflow"
GRAPHS_DIR_POINTSTO="$GRAPH_DIR/pointsto"
GRAPHS_DIR_JAVA_POINTSTO="$GRAPH_DIR/java_pointsto"
GRAMMAR_DIR="./graphs_grammars/grammars"
GRAMMAR_FILE_DATAFLOW="$GRAMMAR_DIR/rules_dataflow.txt"
GRAMMAR_FILE_POINTSTO="$GRAMMAR_DIR/rules_pointsto.txt"
GRAMMAR_FILE_JAVA_POINTSTO="$GRAMMAR_DIR/rules_java_pointsto.txt"

LOG_DIR="./logs/"  # Directory to store logs

# Create the logs directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Get the current date and time for the log file name
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOG_FILE="$LOG_DIR/results_paper_executables_$TIMESTAMP.log"  # Log file with date and time

# Clear or create the log file
echo "Log file for execution runs mentioned in the paper" > "$LOG_FILE"
echo "Run started at $(date)" >> "$LOG_FILE"
echo "***********************************************************" >> "$LOG_FILE"
echo "***********************************************************" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

# List of executables to run (add the names of the executables you want to run from the BIN_DIR)
EXECUTABLES_TO_RUN=( 
 "e-centric-bi-grammar_driven-parallel"
 "v-centric-fw-grammar_driven-label_indexed-parallel"
 "v-centric-bi-grammar_driven-label_indexed-parallel"
 "e-centric-bi-grammar_driven" 
 "v-centric-bi-grammar_driven-label_indexed" 
 "v-centric-fw-grammar_driven-label_indexed" 
 "v-centric-bw-grammar_driven-label_indexed" 
 "v-centric-bi-grammar_driven-label_indexed-vec_copy"
 "e-centric-bi-topo_driven-rule_lookup" 
 "v-centric-bi-topo_driven-sliding_ptrs-rule_lookup"
 )

# Function to process graphs and grammar
run_executables() {
  local GRAPHS_DIR=$1
  local GRAMMAR_FILE=$2

  if [ ! -d "$GRAPHS_DIR" ]; then
    echo "Error: graphs directory $GRAPHS_DIR not found!" | tee -a "$LOG_FILE"
    return
  fi

  if [ ! -f "$GRAMMAR_FILE" ]; then
    echo "Error: grammar file $GRAMMAR_FILE not found!" | tee -a "$LOG_FILE"
    return
  fi

  for executable_name in "${EXECUTABLES_TO_RUN[@]}"; do
    executable="$BIN_DIR/$executable_name"

    if [ -x "$executable" ]; then  # Check if the file is executable
      for graph in "$GRAPHS_DIR"/*; do
        echo "-----------------------------------------------------------"
        echo -e "Executable:\t$executable"
        echo -e "Graph:\t$graph"
        echo -e "Grammar:\t$GRAMMAR_FILE"

        # Run the executable with the graph and grammar file and log the output
        echo -e "Running command:\t$executable $graph $GRAMMAR_FILE" | tee -a "$LOG_FILE"
        "$executable" "$graph" "$GRAMMAR_FILE" >> "$LOG_FILE" 2>&1

        # Capture exit status
        if [ $? -ne 0 ]; then
          echo "Error: $executable failed with graph $graph and grammar $GRAMMAR_FILE." | tee -a "$LOG_FILE"
        else
          echo "$executable ran successfully with graph $graph and grammar $GRAMMAR_FILE."
        fi

         echo "-----------------------------------------------------------"
        echo "" | tee -a "$LOG_FILE"
      done
    else
      echo "Skipping $executable, not an executable file." | tee -a "$LOG_FILE"
    fi
  done
}

echo "-----------------DATAFLOW ANALYSIS STARTED-----------------" | tee -a "$LOG_FILE"
# Execute for dataflow graphs
run_executables "$GRAPHS_DIR_DATAFLOW" "$GRAMMAR_FILE_DATAFLOW"
echo "-----------------DATAFLOW ANALYSIS FINISHED----------------" | tee -a "$LOG_FILE"

echo "-----------------POINTSTO ANALYSIS STARTED-----------------" | tee -a "$LOG_FILE"
# Execute for pointsto graphs
run_executables "$GRAPHS_DIR_POINTSTO" "$GRAMMAR_FILE_POINTSTO"
echo "-----------------POINTSTO ANALYSIS FINISHED-----------------" | tee -a "$LOG_FILE"

echo "-----------------JAVA POINTSTO ANALYSIS STARTED-------------" | tee -a "$LOG_FILE"
# Execute for java_pointsto graphs
run_executables "$GRAPHS_DIR_JAVA_POINTSTO" "$GRAMMAR_FILE_JAVA_POINTSTO"
echo "-----------------JAVA POINTSTO ANALYSIS FINISHED------------" | tee -a "$LOG_FILE"

# Finalize the log
echo "Run completed at $(date)" >> "$LOG_FILE"
echo "All executions finished." | tee -a "$LOG_FILE"
