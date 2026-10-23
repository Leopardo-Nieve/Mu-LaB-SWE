#!/usr/bin/env bash
#SBATCH --time=01:00:00

PROJECT_NAME="convergenceseabed" # maximum 17? characters without underscores? (generates error with write(fpath,...) statement in write_csv

PARAM_FILE="./doc/parameters.cfg"
DATA_DIR="./data"
RESULTS_DIR="./results"
CONV_RESULTS_DIR="${RESULTS_DIR}/convergence"
REFINEMENTS=("$@")
DX_REF=2.0
DT_REF=0.2

# create directory for convergence data
mkdir -p "${CONV_RESULTS_DIR}"

# clean potentially existing files in directory for convergence data
rm -v "${CONV_RESULTS_DIR}"/*

# compile
echo "Compiling executable..."
make
echo "Executable compiled."

# replace project name
	sed -i "s/project_name = [a-zA-Z0-9]*/project_name = ${PROJECT_NAME}/g" "${PARAM_FILE}"


for r in "${REFINEMENTS[@]}"; do
	# update values of dx and dt for current refinement
	DX=$(python -c "print(f'{${DX_REF} / ${r}:.8f}')")
	DT=$(python -c "print(f'{${DT_REF} / (${r}*${r}):.10f}')")

	sed -i "s/dx = [0-9]*\.[0-9]*/dx = ${DX}/g" "${PARAM_FILE}"
	sed -i "s/dt = [0-9]*\.[0-9]*d\?[0-9]*/dt = ${DT}d0/g" "${PARAM_FILE}"
	
	# execute program
	echo "Executing program for refinement r=${r} (dx = ${DX})."
	./bin/main
	echo "Execution complete."

	# locate latest results
	LATEST_DIR=$(ls -t "${DATA_DIR}" 2>/dev/null | head -1)
	if [ -z "${LATEST_DIR}" ]; then
		echo "  ERROR: no directory found in ${DATA_DIR}/"
		exit 1
	fi
	LATEST_DIR=${DATA_DIR}/${LATEST_DIR}

	LATEST_CSV=$(ls -t "${LATEST_DIR}"/* 2>/dev/null | head -1)
	if [ -z "${LATEST_CSV}" ]; then
		echo "  ERROR: no CSV found in ${LATEST_DIR}/"
		exit 1
	fi
	
	# descriptive name to work with python plotting code
	DEST_CSV="${CONV_RESULTS_DIR}/r${r}_dx${DX}.csv"
	
	# move results to convergence analysis
	mv "${LATEST_CSV}" "${DEST_CSV}"

	#debug
	echo ls ./results/convergence:
	ls ./results/convergence

	sleep 1 # ensures new directory name for quick simulations

done
# delete unnecessary files
echo "Cleaning repository..."
make clean
echo "Repository cleaned."

# plot convergence analysis
python3 bin/plot_convergence.py

