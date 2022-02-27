# +--------------------------------------------------------------------------+ #
# Root of mol-infer
# Do not add "/" or "\" to the end of the path (especially on windows).
# Examples:
# Windows    "G:/mol-infer"
# Linux      "~/mol-infer"
# +--------------------------------------------------------------------------+ #
MOLINFER_ROOT="D:/lab/mol-infer"

# +--------------------------------------------------------------------------+ #
# Path to CPLEX executable
# Examples:
# Windows    "G:/CPLEX/CPLEX_Studio1210/cplex/bin/x64_win64/cplex.exe"
# Linux      "/opt/ibm/cplex_12.10/cplex/bin/x86-64_linux/cplex"
# +--------------------------------------------------------------------------+ #
CPLEX_PATH="D:/lab/CPLEX_Studio1210/cplex/bin/x64_win64/cplex.exe"

# +--------------------------------------------------------------------------+ #
# Operating system
# Should be windows / linux / macos
# +--------------------------------------------------------------------------+ #
OS="windows"

# +--------------------------------------------------------------------------+ #
# Default values
# These values are used when you input nothing and just press enter.
# By default, they are set to examples and can be used to demonstrate
# how this project (mol-infer) works.
# In case you have your own prefered parameters, change them accordingly.
# Use ${VAR} to use previously defined variables.
# +--------------------------------------------------------------------------+ #

# Acyclic train
ACYCLIC_DEFAULT_MOLECULES_FILE="${MOLINFER_ROOT}/examples/Acyclic/1Rt_eliminated.sdf"
ACYCLIC_DEFAULT_TARGET_VALUES_FILE="${MOLINFER_ROOT}/examples/Acyclic/1Rt_valuelist.txt"
ACYCLIC_DEFAULT_ANN_LAYERS="10"
ACYCLIC_DEFAULT_TASK_PREFIX="1Rt"

# Acyclic infer
ACYCLIC_DEFAULT_TARGET_VALUE=1900
ACYCLIC_DEFAULT_PARAM_N=15
ACYCLIC_DEFAULT_PARAM_DIA=10
ACYCLIC_DEFAULT_PARAM_K=2
ACYCLIC_DEFAULT_PARAM_DMAX=3
ACYCLIC_DEFAULT_PARAM_BLK=3
ACYCLIC_DEFAULT_PARAM_BHK=2
ACYCLIC_DEFAULT_MODULE_4_A=3600
ACYCLIC_DEFAULT_MODULE_4_B=10000000
ACYCLIC_DEFAULT_MODULE_4_C=100

# Cyclic train
CYCLIC_DEFAULT_MOLECULES_FILE="${MOLINFER_ROOT}/examples/Cyclic/BP/BP.sdf"
CYCLIC_DEFAULT_TARGET_VALUES_FILE="${MOLINFER_ROOT}/examples/Cyclic/BP/BP_values.txt"
CYCLIC_DEFAULT_ANN_LAYERS="20 10"
CYCLIC_DEFAULT_TASK_PREFIX="BP"

# Cyclic infer
CYCLIC_DEFAULT_TARGET_VALUE=100
CYCLIC_DEFAULT_SPEC_FILE="${MOLINFER_ROOT}/examples/Cyclic/chemical_specification/instance_a.txt"
CYCLIC_DEFAULT_MODULE_4_A=10
CYCLIC_DEFAULT_MODULE_4_B=10000000
CYCLIC_DEFAULT_MODULE_4_C=5
CYCLIC_DEFAULT_MODULE_4_D=2

# Cyclic_improved train
CYCLIC_IMPROVED_DEFAULT_MOLECULES_FILE="${MOLINFER_ROOT}/examples/Cyclic_improved/BP/BP.sdf"
CYCLIC_IMPROVED_DEFAULT_TARGET_VALUES_FILE="${MOLINFER_ROOT}/examples/Cyclic_improved/BP/BP_value.csv"
CYCLIC_IMPROVED_DEFAULT_ANN_LAYERS="20 10"
CYCLIC_IMPROVED_DEFAULT_TASK_PREFIX="BP"

# Cyclic_improved infer
CYCLIC_IMPROVED_DEFAULT_TARGET_VALUE=300
CYCLIC_IMPROVED_DEFAULT_SPEC_FILE="${MOLINFER_ROOT}/examples/Cyclic_improved/chemical_specification/instance_a.txt"
CYCLIC_IMPROVED_DEFAULT_MODULE_4_A=10
CYCLIC_IMPROVED_DEFAULT_MODULE_4_B=10000000
CYCLIC_IMPROVED_DEFAULT_MODULE_4_C=5
CYCLIC_IMPROVED_DEFAULT_MODULE_4_D=2

# 2L-model train
TWOLMODEL_DEFAULT_MOLECULES_FILE="${MOLINFER_ROOT}/examples/2L-model/Kow/Kow_eli_C_O_N.sdf"
TWOLMODEL_DEFAULT_TARGET_VALUES_FILE="${MOLINFER_ROOT}/examples/2L-model/Kow/Kow_values.txt"
TWOLMODEL_DEFAULT_ANN_ITERS=10000
TWOLMODEL_DEFAULT_ANN_LAYERS="77"
TWOLMODEL_DEFAULT_TASK_PREFIX="Kow"

# 2L-model infer
TWOLMODEL_DEFAULT_TARGET_VALUE=3.2
TWOLMODEL_DEFAULT_SPEC_FILE="${MOLINFER_ROOT}/examples/2L-model/topological_description/instance_a.txt"
TWOLMODEL_DEFAULT_FRINGE_FILE="${MOLINFER_ROOT}/examples/2L-model/fringe_set/ins_a_fringe.txt"
TWOLMODEL_DEFAULT_MODULE_4_A=10
TWOLMODEL_DEFAULT_MODULE_4_B=10000
TWOLMODEL_DEFAULT_MODULE_4_C=5
TWOLMODEL_DEFAULT_MODULE_4_D=10
TWOLMODEL_DEFAULT_MODULE_4_E=10000
TWOLMODEL_DEFAULT_MODULE_4_F=2

# 2LMM-LLR train
TWOLMMLLR_DEFAULT_MOLECULES_FILE="${MOLINFER_ROOT}/examples/2LMM-LLR/Hc_Hdef.sdf"
TWOLMMLLR_DEFAULT_ATOM_LIMIT=""
TWOLMMLLR_DEFAULT_TARGET_VALUES_FILE="${MOLINFER_ROOT}/examples/2LMM-LLR/Hc_values.txt"
TWOLMMLLR_DEFAULT_LLR_LAMBDA=1.9e-4
TWOLMMLLR_DEFAULT_TASK_PREFIX="Hc"

# 2LMM-LLR infer
TWOLMMLLR_DEFAULT_TARGET_LOWER=1900
TWOLMMLLR_DEFAULT_TARGET_UPPER=1920
TWOLMMLLR_DEFAULT_SPEC_FILE="${MOLINFER_ROOT}/examples/2LMM-LLR/instance_b4_test_2LMM.txt"
TWOLMMLLR_DEFAULT_FRINGE_FILE="${MOLINFER_ROOT}/examples/2LMM-LLR/ins_b4_test_fringe_2LMM.txt"
TWOLMMLLR_DEFAULT_MODULE_4_A=10
TWOLMMLLR_DEFAULT_MODULE_4_B=10000
TWOLMMLLR_DEFAULT_MODULE_4_C=5
TWOLMMLLR_DEFAULT_MODULE_4_D=10
TWOLMMLLR_DEFAULT_MODULE_4_E=10000
TWOLMMLLR_DEFAULT_MODULE_4_F=2

# OS-specific configurations
if [ $OS = "Windows" ] || [ $OS = "windows" ]; then
    UNIX_TOOLS_PATH="${MOLINFER_ROOT}/cygwin/bin/"
    OS="windows"
    PYTHON="${MOLINFER_ROOT}/python-venv/Scripts/python"
    EXE_POSTFIX=".exe"
elif [ $OS = "Linux" ] || [ $OS = "linux" ]; then
    UNIX_TOOLS_PATH=""
    OS="linux"
    PYTHON="${MOLINFER_ROOT}/python-venv/bin/python"
    EXE_POSTFIX=""
elif [ $OS = "MacOS" ] || [ $OS = "macos" ]; then
    UNIX_TOOLS_PATH=""
    OS="macos"
    PYTHON="${MOLINFER_ROOT}/python-venv/bin/python"
    EXE_POSTFIX=""
fi
