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

# First, input prefix for files generated by procedure 'train'.
# Prefix is stored in TASK_PREFIX.
# TEMPLATE_DEFAULT_TASK_PREFIX is used as default prefix. This value is defined
# in mol-infer_config.sh.
# Change default prefix and example as needed.
echo "Please enter the prefix of input files."
echo "E.g.: You have Hc_linreg.txt, then the prefix should be Hc"
echo "Default: ${TEMPLATE_DEFAULT_TASK_PREFIX}"
read -e -p "$(echo -e "[${Green}prefix${Color_Off}]: ")" TASK_PREFIX
TASK_PREFIX=${TASK_PREFIX:-$TEMPLATE_DEFAULT_TASK_PREFIX}
echo ""

# Then, guide users to enter other files and parameters.
# The following 5 lines asks users to input a target value,
# and store this value to TARGET_VALUE.
# TEMPLATE_DEFAULT_TARGET_VALUE is also defined in
# mol-infer_config.sh.
# Duplicate the following 5 lines to add more parameters or input files.
echo "Please enter a target value."
echo "Default: ${TEMPLATE_DEFAULT_TARGET_VALUE}"
read -e -p "$(echo -e "[${Green}lower bound${Color_Off}]: ")" TARGET_VALUE
TARGET_VALUE=${TARGET_VALUE:-$TEMPLATE_DEFAULT_TARGET_VALUE}
echo ""

echo "-----------------------------------------------------"
echo "prefix        ${TASK_PREFIX}"
echo "target value  ${TARGET_VALUE}"
echo "-----------------------------------------------------"

while true; do
    read -p "Proceed? [y/n] " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please enter yes (y) or no (n).";;
    esac
done

# The next step is to generate and execute corresponding commands.
echo -e "${Yellow}"
echo "Solving MILP..."
echo ""

# Example for running python scripts.
# Call python using $PYTHON. This variable points to the python executable in
# python virtual environment created by user.
$PYTHON $MOLINFER_ROOT/2LMM-LLR/Module_3/infer_2LMM_LLR.py "$TASK_PREFIX" \
    $TARGET_LOWER $TARGET_UPPER "$SPEC_FILE" "$FRINGE_FILE" \
    "${TASK_PREFIX}_MILP" "$CPLEX_PATH"

# Check return value.
# Throw error and terminate in case it fails to run.
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "MILP is infesable or error occured while solving."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

# Print necessary information.
echo ""
echo "Generating isomers..."
echo ""

# Example for calling binaries.
# Binary files should be located in $MOLINFER_ROOT/<package_name>/bin/.
$MOLINFER_ROOT/2LMM-LLR/bin/generate_isomers \
    "${TASK_PREFIX}_MILP.sdf" \
    $MODULE_4_A $MODULE_4_B $MODULE_4_C $MODULE_4_D $MODULE_4_E $MODULE_4_F \
    "${TASK_PREFIX}_output.sdf" "${TASK_PREFIX}_MILP_partition.txt" \
    "${FRINGE_FILE}"

# Check return value.
if [ "$?" != "0" ]; then
    echo -e "${Red}"
    echo "Error occured while generating isomers."
    echo -e "${Color_Off}"
    read -p "Press enter to continue."
    exit 1
fi

# Informs users that the program has ended successfully.
echo -e "${Green}"
echo "Done."
echo -e "${Color_Off}"

# clear stdin
while read -r -t 0; do read -r; done
read -p "Press enter to continue."

# Show the output files.
echo ""
echo "Result has been saved to:"
echo "${TASK_PREFIX}_output.sdf"
echo "These file contains generated molecules with desired property."
echo "See documents in Cyclic_improved folder for details about output files."
echo ""

# clear stdin
while read -r -t 0; do read -r; done
read -p "Press enter to continue."