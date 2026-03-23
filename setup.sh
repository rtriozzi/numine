#!/usr/bin/env bash

# base repository
REPO_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# cafana base directory
CAFANA_DIR="${REPO_DIR}/cafana"

# list of analysis directories
analysis_dirs=("cc1e0pi" "cc1e0pi0p" "cc1mu0pi" "michels" "ncpi0")

echo "Repo directory:   $REPO_DIR"
echo "CAFANA directory: $CAFANA_DIR"
echo ""

for dir in "${analysis_dirs[@]}"; do

  analysis_path="${CAFANA_DIR}/${dir}"
  run_dir="${analysis_path}/run"

  if [ ! -d "$analysis_path" ]; then
    echo "Warning: ${analysis_path} does not exist. Creating it."
    mkdir -p "$analysis_path"
  fi

  # create run directories
  mkdir -p "$run_dir"
  mkdir -p "$run_dir/plots"
  mkdir -p "$run_dir/debug"

  # setup env var names
  var_name=$(echo "${dir}_run_dir" | tr '[:lower:]' '[:upper:]')
  export ${var_name}="$run_dir"

  # print info
  echo "Created:  $run_dir"
  echo "Exported: ${var_name}=${run_dir}"
  echo ""
done