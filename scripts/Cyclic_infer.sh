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

echo "------"
echo "Parameters for MILP (Module 3)."
echo "------"
echo ""

echo "Please enter the prefix of input files."
echo "E.g.: You have MP_weights.txt and MP_biases.txt,"
echo "      then the prefix should be MP"
echo "Default: ${CYCLIC_DEFAULT_TASK_PREFIX}"
read -e -p "$(echo -e "[${Green}prefix${Color_Off}]: ")" TASK_PREFIX
TASK_PREFIX=${TASK_PREFIX:-$CYCLIC_DEFAULT_TASK_PREFIX}
echo ""

echo "Please enter a target value."
echo "Default: ${CYCLIC_DEFAULT_TARGET_VALUE}"
read -e -p "$(echo -e "[${Green}target value${Color_Off}]: ")" TARGET_VALUE
TARGET_VALUE=${TARGET_VALUE:-$CYCLIC_DEFAULT_TARGET_VALUE}
echo ""

echo "Please supply a chemical specification file (see documents)."
echo "Default: ${CYCLIC_DEFAULT_SPEC_FILE}"
read -e -p "$(echo -e "[${Green}filename${Color_Off}]: ")" SPEC_FILE
SPEC_FILE=${SPEC_FILE:-$CYCLIC_DEFAULT_SPEC_FILE}
echo ""

echo "------"
echo "Parameters for Graph Generation (Module 4)."
echo "------"
echo ""

echo "Please enter a computation time limit (in seconds)."
echo "Default: ${CYCLIC_DEFAULT_MODULE_4_A}"
read -e -p "$(echo -e "[${Green}time limit${Color_Off}]: ")" MODULE_4_A
MODULE_4_A=${MODULE_4_A:-$CYCLIC_DEFAULT_MODULE_4_A}
echo ""

echo "Please enter a upper bound on the number of stored partial"
echo "feature vectors."
echo "Default: ${CYCLIC_DEFAULT_MODULE_4_B}"
read -e -p "$(echo -e "[${Green}value${Color_Off}]: ")" MODULE_4_B
MODULE_4_B=${MODULE_4_B:-$CYCLIC_DEFAULT_MODULE_4_B}
echo ""

echo "Please enter a upper bound on the number of sample graphs"
echo "stored per feature vector."
echo "Default: ${CYCLIC_DEFAULT_MODULE_4_C}"
read -e -p "$(echo -e "[${Green}value${Color_Off}]: ")" MODULE_4_C
MODULE_4_C=${MODULE_4_C:-$CYCLIC_DEFAULT_MODULE_4_C}
echo ""

echo "Please enter a upper bound on the number of output graphs."
echo "Default: ${CYCLIC_DEFAULT_MODULE_4_D}"
read -e -p "$(echo -e "[${Green}value${Color_Off}]: ")" MODULE_4_D
MODULE_4_D=${MODULE_4_D:-$CYCLIC_DEFAULT_MODULE_4_D}
echo ""

echo -e "${Yellow}"
echo "Solving MILP..."
echo ""

$PYTHON $MOLINFER_ROOT/Cyclic/bin/infer_cyclic_graphs_ec_id.py \
    "$TASK_PREFIX" $TARGET_VALUE "$SPEC_FILE" \
    "${TASK_PREFIX}_MILP" 1 "$CPLEX_PATH"
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "MILP is infesable or error occured while solving."
    read -p "Press enter to continue."
    echo -e "${Color_Off}"
    exit 1
fi

echo ""
echo "Generating isomers..."
echo ""

$MOLINFER_ROOT/Cyclic/bin/$OS/generate_partition \
    "${TASK_PREFIX}_MILP.sdf" \
    "${TASK_PREFIX}_MILP_partition.txt"
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while generating partition."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

$MOLINFER_ROOT/Cyclic/bin/$OS/generate_isomers \
    "${TASK_PREFIX}_MILP.sdf" \
    $MODULE_4_A $MODULE_4_B $MODULE_4_C $MODULE_4_D \
    "${TASK_PREFIX}_output.sdf" \
    "${TASK_PREFIX}_MILP_partition.txt"
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while generating isomers."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

echo -e "${Green}"
echo "Done."
echo -e "${Color_Off}"
echo "Result has been saved to:"
echo "${TASK_PREFIX}_output.sdf"
echo ""
read -p "Press enter to continue"