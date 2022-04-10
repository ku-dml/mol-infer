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

echo "Please supply molecules used for training (sdf format)."
echo "Default: ${TWOLMMLLR_DEFAULT_MOLECULES_FILE}"
read -e -p "$(echo -e "[${Green}molecules file${Color_Off}]: ")" MOLECULES_FILE
MOLECULES_FILE=${MOLECULES_FILE:-$TWOLMMLLR_DEFAULT_MOLECULES_FILE}
echo ""

echo "Please provide atom limitations."
echo "Leave blank to bypass this step."
echo "Example: H C O N"
echo "Default: No limitations"
read -e -p "$(echo -e "[${Green}atoms${Color_Off}]: ")" ATOM_LIMIT
ATOM_LIMIT=${ATOM_LIMIT:-$TWOLMMLLR_DEFAULT_ATOM_LIMIT}
echo ""

echo "Please supply target values for the molecules."
echo "Default: ${TWOLMMLLR_DEFAULT_TARGET_VALUES_FILE}"
read -e -p "$(echo -e "[${Green}values file${Color_Off}]: ")" TARGET_VALUES_FILE
TARGET_VALUES_FILE=${TARGET_VALUES_FILE:-$TWOLMMLLR_DEFAULT_TARGET_VALUES_FILE}
echo ""

echo "Please enter the value of λ to be used in LLR."
echo "Default: ${TWOLMMLLR_DEFAULT_LLR_LAMBDA}"
read -e -p "$(echo -e "[${Green}λ${Color_Off}]: ")" LLR_LAMBDA
LLR_LAMBDA=${LLR_LAMBDA:-$TWOLMMLLR_DEFAULT_LLR_LAMBDA}
echo ""

echo "Please enter the prefix for result files."
echo "Default: ${TWOLMMLLR_DEFAULT_TASK_PREFIX}"
read -e -p "$(echo -e "[${Green}prefix${Color_Off}]: ")" TASK_PREFIX
TASK_PREFIX=${TASK_PREFIX:-$TWOLMMLLR_DEFAULT_TASK_PREFIX}
echo ""

echo -e "${Yellow}"
echo "Calculating molecule descriptors..."
echo ""

$PYTHON $MOLINFER_ROOT/2LMM-LLR/src/Module_1/eliminate.py \
    "$MOLECULES_FILE" "${TASK_PREFIX}_elimanated.sdf"
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while eliminating."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

echo ""
if [ -z "$ATOM_LIMIT" ]
then
    echo "Atom limitation not set. Skip."
    FV_INPUT="${TASK_PREFIX}_elimanated.sdf"
    echo "SDF file for calculating descriptors: $FV_INPUT"
else
    $PYTHON $MOLINFER_ROOT/2LMM-LLR/src/Module_1/limit_atoms.py \
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

$MOLINFER_ROOT/2LMM-LLR/bin/FV_2LMM_V018 "${FV_INPUT}" "${TASK_PREFIX}"
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while calculating descriptors."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

echo ""
echo "Training LLR..."
echo ""

$PYTHON $MOLINFER_ROOT/2LMM-LLR/src/Module_2/lasso_eval_linreg.py \
    "${TASK_PREFIX}_desc_norm.csv" "$TARGET_VALUES_FILE" \
    "${TASK_PREFIX}_linreg.txt" $LLR_LAMBDA
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while training."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

${UNIX_TOOLS_PATH}cp "$TARGET_VALUES_FILE" "${TASK_PREFIX}_values.txt"

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