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
echo "Default: ${TWOLMODEL_DEFAULT_MOLECULES_FILE}"
read -e -p "$(echo -e "[${Green}molecules file${Color_Off}]: ")" MOLECULES_FILE
MOLECULES_FILE=${MOLECULES_FILE:-$TWOLMODEL_DEFAULT_MOLECULES_FILE}
echo ""

echo "Please supply target values for the molecules."
echo "Default: ${TWOLMODEL_DEFAULT_TARGET_VALUES_FILE}"
read -e -p "$(echo -e "[${Green}values file${Color_Off}]: ")" TARGET_VALUES_FILE
TARGET_VALUES_FILE=${TARGET_VALUES_FILE:-$TWOLMODEL_DEFAULT_TARGET_VALUES_FILE}
echo ""

echo "Please enter the number of iterations."
echo "Default: ${TWOLMODEL_DEFAULT_ANN_ITERS}"
read -e -p "$(echo -e "[${Green}iterations${Color_Off}]: ")" ANN_ITERS
ANN_ITERS=${ANN_ITERS:-$TWOLMODEL_DEFAULT_ANN_ITERS}
echo ""

echo "Please enter the hidden sizes of ANN."
echo "Use space to seperate numbers."
echo "Default: ${TWOLMODEL_DEFAULT_ANN_LAYERS}"
read -e -p "$(echo -e "[${Green}layer sizes${Color_Off}]: ")" ANN_LAYERS
ANN_LAYERS=${ANN_LAYERS:-$TWOLMODEL_DEFAULT_ANN_LAYERS}
echo ""

echo "Please enter the prefix for result files."
echo "E.g.: If the prefix is Kow, you will get "
echo "      Kow_weights.txt, Kow_biases.txt, etc."
echo "Default: ${TWOLMODEL_DEFAULT_TASK_PREFIX}"
read -e -p "$(echo -e "[${Green}prefix${Color_Off}]: ")" TASK_PREFIX
TASK_PREFIX=${TASK_PREFIX:-$TWOLMODEL_DEFAULT_TASK_PREFIX}
echo ""

echo -e "${Yellow}"
echo "Calculating molecule descriptors..."
echo ""

$PYTHON $MOLINFER_ROOT/2L-model/Module_1/files/eliminate.py \
    "$MOLECULES_FILE" "${TASK_PREFIX}_elimanated.sdf"
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while eliminating."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

$MOLINFER_ROOT/2L-model/bin/2L_FV \
    "${TASK_PREFIX}_elimanated.sdf" "${TASK_PREFIX}"
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

$PYTHON $MOLINFER_ROOT/2L-model/Module_2/files/2L_ANN.py \
    "${TASK_PREFIX}_desc_norm.csv" "$TARGET_VALUES_FILE" \
    "$TASK_PREFIX" $ANN_ITERS $ANN_LAYERS
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
echo "${TASK_PREFIX}_desc_norm.csv"
echo "${TASK_PREFIX}_fringe.txt"
echo "${TASK_PREFIX}_weights.txt"
echo "${TASK_PREFIX}_biases.txt"
echo ""
read -p "Press enter to continue"