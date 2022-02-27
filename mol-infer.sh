# ============================================
# A quickstart shell script for mol-infer
# Discrete Mathematics Lab, Kyoto University
# December 2021
# ============================================

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

print_header () {
    echo "====================================================="
    echo "mol-infer"
    echo "Molecular inferring project"
    echo "By Discrete Mathematics Lab, Kyoto University."
    echo "====================================================="
    echo ""
    echo "-----------------------------------------------------"
    echo "Current configuration:"
    echo "mol-infer root   = $MOLINFER_ROOT"
    echo "CPLEX path       = $CPLEX_PATH"
    echo "Operating System = $OS"
    echo "-----------------------------------------------------"
}

# Read configuration
if [ -f "./mol-infer_config.sh" ]; then
    . ./mol-infer_config.sh
else
    echo "Missing mol-infer_config.sh"
    exit 1
fi

if [ $OS = "Windows" ] || [ $OS = "windows" ]; then
    UNIX_TOOLS_PATH="${MOLINFER_ROOT}/cygwin/bin/"
elif [ $OS = "Linux" ] || [ $OS = "linux" ]; then
    UNIX_TOOLS_PATH=""
else
    echo "Unknown Operating system $OS."
    exit 1
fi

while true; do

    ${UNIX_TOOLS_PATH}clear
    print_header

    echo ""
    echo "* Please read documents in each package for "
    echo "* information about the input files and parameters."
    echo ""

    # Select package
    PACKAGE_NAME=""
    PS3="Please select a package: "
    options=("Acyclic" "Cyclic" "Cyclic_improved" "2L-model" "2LMM-LLR" "Quit")
    select opt in "${options[@]}"
    do
        case $opt in
            "Acyclic")
                PACKAGE_NAME="Acyclic"
                break
                ;;
            "Cyclic")
                PACKAGE_NAME="Cyclic"
                break
                ;;
            "Cyclic_improved")
                PACKAGE_NAME="Cyclic_improved"
                break
                ;;
            "2L-model")
                PACKAGE_NAME="2L-model"
                break
                ;;
            "2LMM-LLR")
                PACKAGE_NAME="2LMM-LLR"
                break
                ;;
            "Quit")
                break 2
                ;;
            *) echo "Invalid option $REPLY";;
        esac
    done

    while true; do

        ${UNIX_TOOLS_PATH}clear
        print_header

        echo ""
        echo "-----------------------------------------------------"
        echo "Package ${PACKAGE_NAME}"
        echo "-----------------------------------------------------"
        echo ""

        # Select procedure
        PS3="Please select a procedure: "
        options=("Train" "Infer" "Back")
        select opt in "${options[@]}"
        do
            case $opt in
                "Train")
                    echo ""
                    echo "-----------------------------------------------------"
                    echo "Train"
                    echo "See Module 1 & 2."
                    echo "-----------------------------------------------------"
                    echo ""
                    ${UNIX_TOOLS_PATH}bash ${MOLINFER_ROOT}/scripts/${PACKAGE_NAME}_train.sh
                    break
                    ;;
                "Infer")
                    echo ""
                    echo "-----------------------------------------------------"
                    echo "Infer"
                    echo "See Module 3 & 4."
                    echo "-----------------------------------------------------"
                    echo ""
                    ${UNIX_TOOLS_PATH}bash ${MOLINFER_ROOT}/scripts/${PACKAGE_NAME}_infer.sh
                    break
                    ;;
                "Back")
                    break 2
                    ;;
                *) echo "Invalid option $REPLY";;
            esac
        done
    done
done
