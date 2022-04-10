#!/bin/bash

# Text Colors
Color_Off='\033[0m'       # Text Reset
Black='\033[0;30m'        # Black
Red='\033[0;31m'          # Red
Green='\033[0;32m'        # Green
Yellow='\033[0;33m'       # Yellow
Blue='\033[0;34m'         # Blue
Purple='\033[0;35m'       # Purple
Cyan='\033[0;36m'         # Cyan
White='\033[0;37m'        # White

# Read configuration
. ./mol-infer_config.sh
. ./mol-infer_default_values.sh

# OS-specific configurations
if [ $OS = "Windows" ] || [ $OS = "windows" ]; then
    OS="windows"
    PYTHON="${MOLINFER_ROOT}/python-venv/Scripts/python"
elif [ $OS = "Linux" ] || [ $OS = "linux" ]; then
    OS="linux"
    PYTHON="${MOLINFER_ROOT}/python-venv/bin/python"
elif [ $OS = "MacOS" ] || [ $OS = "macos" ]; then
    OS="macos"
    PYTHON="${MOLINFER_ROOT}/python-venv/bin/python"
fi

# Solver configurations
#if [ $SOLVER_TYPE = "CPLEX" ]; then
#    SOLVER_INDICATOR=1
#    SOLVER_PARAM=$CPLEX_PATH
#elif [ $SOLVER_TYPE = "NEOS" ]; then
#    SOLVER_INDICATOR=3
#    SOLVER_PARAM=$NEOS_EMAIL_ADDR
#fi

SOLVER_INDICATOR=1
SOLVER_PARAM=$CPLEX_PATH

# Guide users to enter other files and parameters.
# The following 5 lines asks users to input the file containing input molecules
# and store the filename to MOLECULES_FILE.
# TEMPLATE_DEFAULT_MOLECULES_FILE is used as default filename and it is defined
# in experiments/mol-infer_config.sh.
# Duplicate the following 5 lines to add more parameters or input files.
echo "Please supply molecules used for training (sdf format)."
echo "Default: ${TEMPLATE_DEFAULT_MOLECULES_FILE}"
read -e -p "$(echo -e "[${Green}molecules file${Color_Off}]: ")" MOLECULES_FILE
MOLECULES_FILE=${MOLECULES_FILE:-$TEMPLATE_DEFAULT_MOLECULES_FILE}
echo ""

# Ask users to provide a prefix for output files.
# This prefix will be added to all output files and it is also the identifier
# used in the next stage (infer) to search for input files.
# Default value is defined in experiments/mol-infer_config.sh.
echo "Please enter the prefix for result files."
echo "Default: ${TEMPLATE_DEFAULT_TASK_PREFIX}"
read -e -p "$(echo -e "[${Green}prefix${Color_Off}]: ")" TASK_PREFIX
TASK_PREFIX=${TASK_PREFIX:-$TEMPLATE_DEFAULT_TASK_PREFIX}
echo ""

# Generate and execute corresponding commands.
echo -e "${Yellow}"
echo "Calculating molecule descriptors..."
echo ""

# Example for running python scripts.
# Call python using $PYTHON. This variable points to the python executable in
# python virtual environment created by user.
$PYTHON $MOLINFER_ROOT/2LMM-LLR/Module_1/eliminate.py \
    "$MOLECULES_FILE" "${TASK_PREFIX}_elimanated.sdf"

# Check return value.
# Throw error and terminate in case it fails to run.
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while eliminating."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

# Use if statement if some commands are not always executed.
echo ""
if [ -z "$ATOM_LIMIT" ]
then
    echo "Atom limitation not set. Skip."
    FV_INPUT="${TASK_PREFIX}_elimanated.sdf"
    echo "SDF file for calculating descriptors: $FV_INPUT"
else
    $PYTHON $MOLINFER_ROOT/2LMM-LLR/Module_1/limit_atoms.py \
        "${TASK_PREFIX}_elimanated.sdf" $ATOM_LIMIT
    if [ "$?" != "0" ]; then
        echo -e "${Red}"
        echo "Error occured while limiting atoms."
        echo -e "${Color_Off}"
        read -p "Press enter to continue."
        exit 1
    fi
    FV_INPUT="${TASK_PREFIX}_elimanated_${ATOM_LIMIT// /_}.sdf"
    echo "SDF file for calculating descriptors: $FV_INPUT"
fi
echo ""

# Example for calling binaries.
# Binary files corresponding to users' operating system should be in
# $MOLINFER_ROOT/<package_name>/bin/$OS/.
$MOLINFER_ROOT/2LMM-LLR/bin/FV_2LMM_V018 "${FV_INPUT}" "${TASK_PREFIX}"

# Check return value.
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while calculating descriptors."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

# Call unix tools from ${UNIX_TOOLS_PATH}.
${UNIX_TOOLS_PATH}cp "$TARGET_VALUES_FILE" "${TASK_PREFIX}_values.txt"

# Informs users that the program has ended successfully.
# Show the location of output files.
echo -e "${Green}"
echo "Done."
echo -e "${Color_Off}"
echo "Results have been saved to:"
echo "${TASK_PREFIX}_eliminated.sdf"
echo "${TASK_PREFIX}_desc.csv"
echo "${TASK_PREFIX}_desc_norm.csv"
echo "${TASK_PREFIX}_fringe.txt"
echo "${TASK_PREFIX}_linreg.txt"
echo "${TASK_PREFIX}_values.txt (copied from ${TARGET_VALUES_FILE})"
echo ""
read -p "Press enter to continue."