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

echo "------"
echo "Parameters for MILP (Module 3)."
echo "------"
echo ""

echo "Please enter the prefix of input files."
echo "Default: ${TWOLMODEL_DEFAULT_TASK_PREFIX}"
read -e -p "$(echo -e "[${Green}prefix${Color_Off}]: ")" TASK_PREFIX
TASK_PREFIX=${TASK_PREFIX:-$TWOLMODEL_DEFAULT_TASK_PREFIX}
echo ""

echo "Please enter a target value."
echo "Default: ${TWOLMODEL_DEFAULT_TARGET_VALUE}"
read -e -p "$(echo -e "[${Green}target value${Color_Off}]: ")" TARGET_VALUE
TARGET_VALUE=${TARGET_VALUE:-$TWOLMODEL_DEFAULT_TARGET_VALUE}
echo ""

echo "Please supply chemical specification file."
echo "Default: ${TWOLMODEL_DEFAULT_SPEC_FILE}"
read -e -p "$(echo -e "[${Green}filename${Color_Off}]: ")" SPEC_FILE
SPEC_FILE=${SPEC_FILE:-$TWOLMODEL_DEFAULT_SPEC_FILE}
echo ""

echo "Please supply fringe tree file."
echo "Default: ${TWOLMODEL_DEFAULT_FRINGE_FILE}"
read -e -p "$(echo -e "[${Green}filename${Color_Off}]: ")" FRINGE_FILE
FRINGE_FILE=${FRINGE_FILE:-$TWOLMODEL_DEFAULT_FRINGE_FILE}
echo ""

echo "------"
echo "Parameters for Graph Generation (Module 4)."
echo "------"
echo ""

echo "Please enter time limit (in seconds) on each stages of the program."
echo "Default: ${TWOLMODEL_DEFAULT_MODULE_4_A}"
read -e -p "$(echo -e "[${Green}value${Color_Off}]: ")" MODULE_4_A
MODULE_4_A=${MODULE_4_A:-$TWOLMODEL_DEFAULT_MODULE_4_A}
echo ""

echo "Please enter a upper bound on the number"
echo "of stored partial feature vectors."
echo "Default: ${TWOLMODEL_DEFAULT_MODULE_4_B}"
read -e -p "$(echo -e "[${Green}value${Color_Off}]: ")" MODULE_4_B
MODULE_4_B=${MODULE_4_B:-$TWOLMODEL_DEFAULT_MODULE_4_B}
echo ""

echo "Please enter a upper bound on the number"
echo "of graphs stored per base vertex or edge."
echo "Default: ${TWOLMODEL_DEFAULT_MODULE_4_C}"
read -e -p "$(echo -e "[${Green}value${Color_Off}]: ")" MODULE_4_C
MODULE_4_C=${MODULE_4_C:-$TWOLMODEL_DEFAULT_MODULE_4_C}
echo ""

echo "Please enter time limit (in seconds) for enumeration of paths."
echo "Default: ${TWOLMODEL_DEFAULT_MODULE_4_D}"
read -e -p "$(echo -e "[${Green}value${Color_Off}]: ")" MODULE_4_D
MODULE_4_D=${MODULE_4_D:-$TWOLMODEL_DEFAULT_MODULE_4_D}
echo ""

echo "Please enter a upper bound on the number"
echo "of total paths stored during the computation."
echo "Default: ${TWOLMODEL_DEFAULT_MODULE_4_E}"
read -e -p "$(echo -e "[${Green}value${Color_Off}]: ")" MODULE_4_E
MODULE_4_E=${MODULE_4_E:-$TWOLMODEL_DEFAULT_MODULE_4_E}
echo ""

echo "Please enter a upper bound on the number"
echo "of output graphs."
echo "Default: ${TWOLMODEL_DEFAULT_MODULE_4_F}"
read -e -p "$(echo -e "[${Green}value${Color_Off}]: ")" MODULE_4_F
MODULE_4_F=${MODULE_4_F:-$TWOLMODEL_DEFAULT_MODULE_4_F}
echo ""

echo "-----------------------------------------------------"
echo "prefix                        ${TASK_PREFIX}"
echo "target value                  ${TARGET_VALUE}"
echo "chem specification file       ${SPEC_FILE}"
echo "fringe tree file              ${FRINGE_FILE}"
echo "time limit for each stage     ${MODULE_4_A}"
echo "num stored feature vectors    ${MODULE_4_B}"
echo "num graphs per vertex / edge  ${MODULE_4_C}"
echo "time limit for enum paths     ${MODULE_4_D}"
echo "num stored paths              ${MODULE_4_E}"
echo "num output graphs             ${MODULE_4_F}"
echo "-----------------------------------------------------"

while true; do
    read -p "Proceed? [y/n] " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please enter yes (y) or no (n).";;
    esac
done

echo -e "${Yellow}"
echo "Solving MILP..."
echo ""

$PYTHON $MOLINFER_ROOT/2L-model/Module_3/files/infer_graph_2L_fc.py \
    "$TASK_PREFIX" $TARGET_VALUE "$SPEC_FILE" "$FRINGE_FILE" \
    "${TASK_PREFIX}_MILP" 1 "$CPLEX_PATH"
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "MILP is infesable or error occured while solving."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

echo ""
echo "Generating isomers..."
echo ""

$MOLINFER_ROOT/2L-model/bin/generate_isomers \
    "${TASK_PREFIX}_MILP.sdf" \
    $MODULE_4_A $MODULE_4_B $MODULE_4_C $MODULE_4_D $MODULE_4_E $MODULE_4_F \
    "${TASK_PREFIX}_output.sdf" \
    "${TASK_PREFIX}_MILP_partition.txt" \
    "$FRINGE_FILE"
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

# clear stdin
while read -r -t 0; do read -r; done
read -p "Press enter to continue."

echo ""
echo "Result has been saved to:"
echo "${TASK_PREFIX}_output.sdf"
echo "These file contains generated molecules with desired property."
echo "See documents in 2L-model folder for details about output files."
echo ""

# clear stdin
while read -r -t 0; do read -r; done
read -p "Press enter to continue."