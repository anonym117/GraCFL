#!/bin/bash

# Set the paths for the executables and graphs directories, and grammar file
BIN_DIR="./build/bin"
GRAPH_DIR="$HOME/graspan-research/Graphs/ApacheHttpd2.2.18Dataflow/"
GRAMMAR_FILE="$HOME/graspan-research/GrammarFiles/rules_dataflow"  # The first argument is the grammar file

LOG_DIR="./results-log"  # Directory to store logs

# Create the logs directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Get the current date and time for the log file name
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOG_FILE="$LOG_DIR/run_results_log_$TIMESTAMP.txt"  # Log file with date and time

# Clear or create the log file
echo "Log file for execution runs" > "$LOG_FILE"
echo "Run started at $(date)" >> "$LOG_FILE"
echo "----------------------------------------" >> "$LOG_FILE"
echo "----------------------------------------" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

# Check if directories and grammar file exist
if [ ! -d "$BIN_DIR" ]; then
  echo "Error: bin directory not found!" | tee -a "$LOG_FILE"
  exit 1
fi

if [ ! -d "$GRAPH_DIR" ]; then
  echo "Error: graphs directory not found!" | tee -a "$LOG_FILE"
  exit 1
fi

if [ ! -f "$GRAMMAR_FILE" ]; then
  echo "Error: grammar file $GRAMMAR_FILE not found!" | tee -a "$LOG_FILE"
  exit 1
fi

# Loop through all executables in the bin directory
for executable in "$BIN_DIR"/*; do
  if [ -x "$executable" ]; then  # Check if the file is executable
    # Print an empty line
    #echo "" | tee -a "$LOG_FILE"

    echo -e "Executable:\t$executable" | tee -a "$LOG_FILE"

    # Loop through all graph files in the graphs directory
    for graph in "$GRAPH_DIR"/*; do
      echo -e "Graph:\t$graph" 
      echo -e "Grammar:\t$GRAMMAR_FILE" 

      # Run the executable with the graph and grammar file and log the output
      echo -e "Running command:\t$executable $graph $GRAMMAR_FILE" | tee -a "$LOG_FILE"
      "$executable" "$graph" "$GRAMMAR_FILE" >> "$LOG_FILE" 2>&1

      # Capture exit status
      if [ $? -ne 0 ]; then
        echo "Error: $executable failed with graph $graph and grammar $GRAMMAR_FILE." | tee -a "$LOG_FILE"
      else
        echo "$executable ran successfully with graph $graph and grammar $GRAMMAR_FILE." | tee -a "$LOG_FILE"
      fi

      echo "----------------------------------------" | tee -a "$LOG_FILE"

      # Print an empty line
      echo "" | tee -a "$LOG_FILE"
    done
  else
    echo "Skipping $executable, not an executable file." | tee -a "$LOG_FILE"
  fi
done

# Finalize the log
echo "Run completed at $(date)" >> "$LOG_FILE"
echo "All executions finished." | tee -a "$LOG_FILE"
