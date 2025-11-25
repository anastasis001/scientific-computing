# Sample simulator to Dakota system call script

# The first and second command line arguments to the script are the
# names of the Dakota parameters and results files.
params=$1
results=$2

# ---------
# EXECUTION
# ---------

python3 wrappercurrent.py
