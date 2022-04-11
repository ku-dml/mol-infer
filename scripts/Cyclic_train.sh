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
elif [ $OS = "MacOS" ] || [ $OS = "macOS" ] || [ $OS = "macos" ]; then
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

echo "Please supply molecules used for training (sdf format)."
echo "Default: ${CYCLIC_DEFAULT_MOLECULES_FILE}"
read -e -p "$(echo -e "[${Green}molecules file${Color_Off}]: ")" MOLECULES_FILE
MOLECULES_FILE=${MOLECULES_FILE:-$CYCLIC_DEFAULT_MOLECULES_FILE}
echo ""

echo "Please supply target values for the molecules."
echo "Default: ${CYCLIC_DEFAULT_TARGET_VALUES_FILE}"
read -e -p "$(echo -e "[${Green}values file${Color_Off}]: ")" TARGET_VALUES_FILE
TARGET_VALUES_FILE=${TARGET_VALUES_FILE:-$CYCLIC_DEFAULT_TARGET_VALUES_FILE}
echo ""

echo "Please enter the hidden sizes of ANN."
echo "Use space to seperate numbers."
echo "Default: ${CYCLIC_DEFAULT_ANN_LAYERS}"
read -e -p "$(echo -e "[${Green}layer sizes${Color_Off}]: ")" ANN_LAYERS
ANN_LAYERS=${ANN_LAYERS:-$CYCLIC_DEFAULT_ANN_LAYERS}
echo ""

echo "Please enter the prefix for result files."
echo "E.g.: If the prefix is MP,"
echo "      you will get MP_weights.txt and MP_biases.txt."
echo "Default: ${CYCLIC_DEFAULT_TASK_PREFIX}"
read -e -p "$(echo -e "[${Green}prefix${Color_Off}]: ")" TASK_PREFIX
TASK_PREFIX=${TASK_PREFIX:-$CYCLIC_DEFAULT_TASK_PREFIX}
echo ""

echo -e "${Yellow}"
echo "Calculating molecule descriptors..."
echo ""

$MOLINFER_ROOT/Cyclic/bin/CHECKER "$MOLECULES_FILE"
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Cyclic check failed."
    echo "Please check the document for details."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

$PYTHON $MOLINFER_ROOT/Cyclic/src/Module_1/eliminate.py \
    "$MOLECULES_FILE" "${TASK_PREFIX}_elimanated.sdf"
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while eliminating."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

$MOLINFER_ROOT/Cyclic/bin/FV_ec \
    "${TASK_PREFIX}_elimanated.sdf" "${TASK_PREFIX}_desc.csv"
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while calculating descriptors."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

echo ""
echo "Training ANN..."
echo ""

$PYTHON $MOLINFER_ROOT/Cyclic/src/Module_2/scikit_chemgraph_learning_lim.py \
    "${TASK_PREFIX}_desc.csv" "$TARGET_VALUES_FILE" \
    "$TASK_PREFIX" $ANN_LAYERS
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while training."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

echo -e "${Green}"
echo "Done."
echo -e "${Color_Off}"
echo "Results have been saved to:"
echo "${TASK_PREFIX}_eliminated.sdf"
echo "${TASK_PREFIX}_desc.csv"
echo "${TASK_PREFIX}_weights.txt"
echo "${TASK_PREFIX}_biases.txt"
echo ""
read -p "Press enter to continue."