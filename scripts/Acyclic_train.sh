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

echo "Please supply molecules used for training (sdf format)."
echo "Default: ${ACYCLIC_DEFAULT_MOLECULES_FILE}"
read -e -p "$(echo -e "[${Green}molecules file${Color_Off}]: ")" MOLECULES_FILE
MOLECULES_FILE=${MOLECULES_FILE:-$ACYCLIC_DEFAULT_MOLECULES_FILE}
echo ""

echo "Please supply target values for the molecules."
echo "Default: ${ACYCLIC_DEFAULT_TARGET_VALUES_FILE}"
read -e -p "$(echo -e "[${Green}values file${Color_Off}]: ")" TARGET_VALUES_FILE
TARGET_VALUES_FILE=${TARGET_VALUES_FILE:-$ACYCLIC_DEFAULT_TARGET_VALUES_FILE}
echo ""

echo "Please enter the hidden sizes of ANN."
echo "Use space to seperate numbers."
echo "Default: ${ACYCLIC_DEFAULT_ANN_LAYERS}"
read -e -p "$(echo -e "[${Green}layer sizes${Color_Off}]: ")" ANN_LAYERS
ANN_LAYERS=${ANN_LAYERS:-$ACYCLIC_DEFAULT_ANN_LAYERS}
echo ""

echo "Please enter the prefix of result files."
echo "E.g.: If the prefix is 1Rt,"
echo "      you will get 1Rt_weights.txt and 1Rt_biases.txt."
echo "Default: ${ACYCLIC_DEFAULT_TASK_PREFIX}"
read -e -p "$(echo -e "[${Green}prefix${Color_Off}]: ")" TASK_PREFIX
TASK_PREFIX=${TASK_PREFIX:-$ACYCLIC_DEFAULT_TASK_PREFIX}
echo ""

echo -e "${Yellow}"
echo "Calculating molecule descriptors..."
echo ""

$MOLINFER_ROOT/Acyclic/bin/$OS/fv4_in_ex \
    "$MOLECULES_FILE" "${TASK_PREFIX}_desc.csv"
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

$PYTHON $MOLINFER_ROOT/Acyclic/bin/scikit_chemgraph_learning.py \
    "${TASK_PREFIX}_desc.csv" "$TARGET_VALUES_FILE" \
    "$TASK_PREFIX" $ANN_LAYERS
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while training."
    echo -e "${Color_Off}"
    read -p "Press enter to continue"
    exit 1
fi

echo -e "${Green}"
echo "Done."
echo -e "${Color_Off}"
echo "Results have been saved to:"
echo "${TASK_PREFIX}_weights.txt"
echo "${TASK_PREFIX}_biases.txt"
echo "${TASK_PREFIX}_desc.csv"
echo ""
read -p "Press enter to continue."