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

echo "------"
echo "Parameters for MILP (Module 3)."
echo "------"
echo ""

echo "Please enter the prefix of input files."
echo "E.g.: You have 1Rt_weights.txt and 1Rt_biases.txt,"
echo "      then the prefix should be 1Rt"
echo "Default: ${ACYCLIC_DEFAULT_TASK_PREFIX}"
read -e -p "$(echo -e "[${Green}prefix${Color_Off}]: ")" TASK_PREFIX
TASK_PREFIX=${TASK_PREFIX:-$ACYCLIC_DEFAULT_TASK_PREFIX}
echo ""

echo "Please enter a target value."
echo "Default: ${ACYCLIC_DEFAULT_TARGET_VALUE}"
read -e -p "$(echo -e "[${Green}target value${Color_Off}]: ")" TARGET_VALUE
TARGET_VALUE=${TARGET_VALUE:-$ACYCLIC_DEFAULT_TARGET_VALUE}
echo ""

echo "Please enter the number of vertices n*."
echo "Default: ${ACYCLIC_DEFAULT_PARAM_N}"
read -e -p "$(echo -e "[${Green}n*${Color_Off}]: ")" PARAM_N
PARAM_N=${PARAM_N:-$ACYCLIC_DEFAULT_PARAM_N}
echo ""

echo "Please enter the diameter dia*."
echo "Default: ${ACYCLIC_DEFAULT_PARAM_DIA}"
read -e -p "$(echo -e "[${Green}dia*${Color_Off}]: ")" PARAM_DIA
PARAM_DIA=${PARAM_DIA:-$ACYCLIC_DEFAULT_PARAM_DIA}
echo ""

echo "Please enter branch-height parameter k*."
echo "Default: ${ACYCLIC_DEFAULT_PARAM_K}"
read -e -p "$(echo -e "[${Green}k*${Color_Off}]: ")" PARAM_K
PARAM_K=${PARAM_K:-$ACYCLIC_DEFAULT_PARAM_K}
echo ""

echo "Please enter maximum degree d_max of a vertex in the target graph."
echo "This value should be either 3 or 4."
echo "Default: ${ACYCLIC_DEFAULT_PARAM_DMAX}"
read -e -p "$(echo -e "[${Green}d_max${Color_Off}]: ")" PARAM_DMAX
PARAM_DMAX=${PARAM_DMAX:-$ACYCLIC_DEFAULT_PARAM_DMAX}
echo ""

echo "Please enter the k*-branch-leaf-number bl_k*."
echo "Default: ${ACYCLIC_DEFAULT_PARAM_BLK}"
read -e -p "$(echo -e "[${Green}bl_k*${Color_Off}]: ")" PARAM_BLK
PARAM_BLK=${PARAM_BLK:-$ACYCLIC_DEFAULT_PARAM_BLK}
echo ""

echo "Please enter the k*-branch-height bh_k*."
echo "Default: ${ACYCLIC_DEFAULT_PARAM_BHK}"
read -e -p "$(echo -e "[${Green}bh_k*${Color_Off}]: ")" PARAM_BHK
PARAM_BHK=${PARAM_BHK:-$ACYCLIC_DEFAULT_PARAM_BHK}
echo ""

echo "------"
echo "Parameters for Graph Generation (Module 4)."
echo "------"
echo ""

echo "Please enter a computation time limit (in seconds)."
echo "Default: ${ACYCLIC_DEFAULT_MODULE_4_A}"
read -e -p "$(echo -e "[${Green}time limit${Color_Off}]: ")" MODULE_4_A
MODULE_4_A=${MODULE_4_A:-$ACYCLIC_DEFAULT_MODULE_4_A}
echo ""

echo "Please enter a upper bound on the number of stored partial"
echo "feature vectors."
echo "Default: ${ACYCLIC_DEFAULT_MODULE_4_B}"
read -e -p "$(echo -e "[${Green}number features${Color_Off}]: ")" MODULE_4_B
MODULE_4_B=${MODULE_4_B:-$ACYCLIC_DEFAULT_MODULE_4_B}
echo ""

echo "Please enter a upper bound on the number of output graphs."
echo "Default: ${ACYCLIC_DEFAULT_MODULE_4_C}"
read -e -p "$(echo -e "[${Green}number graphs${Color_Off}]: ")" MODULE_4_C
MODULE_4_C=${MODULE_4_C:-$ACYCLIC_DEFAULT_MODULE_4_C}
echo ""

echo -e "${Yellow}"
echo "Solving MILP..."
echo ""

$PYTHON $MOLINFER_ROOT/Acyclic/Module_3/infer_acyclic_graphs.py $TARGET_VALUE \
    $PARAM_N $PARAM_DIA $PARAM_K \
    $PARAM_DMAX $PARAM_BLK $PARAM_BHK \
    $SOLVER_INDICATOR $TASK_PREFIX "$SOLVER_PARAM" "${TASK_PREFIX}_MILP"
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "MILP is infesable or error occured while solving."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

echo ""
echo "2 branches..."
echo ""

$MOLINFER_ROOT/Acyclic/bin/2-branches \
    "${TASK_PREFIX}_MILP.txt" $MODULE_4_A $MODULE_4_B $MODULE_4_C \
    "${TASK_PREFIX}_output_2_branches.sdf" "${TASK_PREFIX}_MILP.sdf"
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while generating isomers."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

echo ""
echo "3 branches..."
echo ""

$MOLINFER_ROOT/Acyclic/bin/3-branches \
    "${TASK_PREFIX}_MILP.txt" $MODULE_4_A $MODULE_4_B $MODULE_4_C \
    "${TASK_PREFIX}_output_3_branches.sdf" "${TASK_PREFIX}_MILP.sdf"
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
echo "${TASK_PREFIX}_output_2_branches.sdf"
echo "${TASK_PREFIX}_output_3_branches.sdf"
echo ""
read -p "Press enter to continue."