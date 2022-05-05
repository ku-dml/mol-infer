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

# Read arguments
MOLINFER_ROOT="$1"
PYTHON="$2"

# Read configuration
. ${MOLINFER_ROOT}/mol-infer_config.sh
. ${MOLINFER_ROOT}/mol-infer_default_values.sh

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
# in mol-infer_config.sh.
# Duplicate the following 5 lines to add more parameters or input files.
echo "Please supply molecules used for training (sdf format)."
echo "Default: ${TEMPLATE_DEFAULT_MOLECULES_FILE}"
read -e -p "$(echo -e "[${Green}molecules file${Color_Off}]: ")" MOLECULES_FILE
MOLECULES_FILE=${MOLECULES_FILE:-$TEMPLATE_DEFAULT_MOLECULES_FILE}
echo ""

# Ask users to provide a prefix for output files.
# This prefix will be added to all output files and it is also the identifier
# used in the next stage (infer) to search for input files.
# Default value is defined in mol-infer_config.sh.
echo "Please enter the prefix for result files."
echo "Default: ${TEMPLATE_DEFAULT_TASK_PREFIX}"
read -e -p "$(echo -e "[${Green}prefix${Color_Off}]: ")" TASK_PREFIX
TASK_PREFIX=${TASK_PREFIX:-$TEMPLATE_DEFAULT_TASK_PREFIX}
echo ""

# The following lines prints out all values for users to double check
echo "-----------------------------------------------------"
echo "molecules file  ${MOLECULES_FILE}"
echo "prefix          ${TASK_PREFIX}"
echo "-----------------------------------------------------"

while true; do
    read -p "Proceed? [y/n] " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please enter yes (y) or no (n).";;
    esac
done

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
# Binary files should be located in $MOLINFER_ROOT/<package_name>/bin/.
$MOLINFER_ROOT/2LMM-LLR/bin/FV_2LMM_V018 "${FV_INPUT}" "${TASK_PREFIX}"

# Check return value.
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while calculating descriptors."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

# Additional operations
cp "$TARGET_VALUES_FILE" "${TASK_PREFIX}_values.txt"

# Informs users that the program has ended successfully.
echo -e "${Green}"
echo "Done."
echo -e "${Color_Off}"

# clear stdin
while read -r -t 0; do read -r; done
read -p "Press enter to continue."

# Show the output files.
echo ""
echo "Results have been saved to:"
echo "${TASK_PREFIX}_eliminated.sdf  This file contains molecules that satisfy condition (ii) in documents."
echo "${TASK_PREFIX}_desc.csv        This file contains descriptors for input molecules."
echo "${TASK_PREFIX}_weights.txt     This file contains weights for the trained ANN."
echo "${TASK_PREFIX}_biases.txt      This file contains biases for the trained ANN."
echo "These files will be used in the [Infer] part."
echo "See documents in Cylic folder for details about output files."
echo ""

# clear stdin
while read -r -t 0; do read -r; done
read -p "Press enter to continue."