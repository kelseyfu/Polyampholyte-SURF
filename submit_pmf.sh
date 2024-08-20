source "inputs/double_chain_0.8_0_0.sh"

# Iterate over lower bounds 0 to 25 in increments of 5
for LBOUND in $(seq 0 5 25)
do
    cp run_pmf.sh run_pmf_$LBOUND.sh
    sed -i "s/LBOUND/$LBOUND/g" run_pmf_$LBOUND.sh
 

    # Define the upper bound
    export UBOUND=$(($LBOUND + 5))
    export LBOUND
    # Run the simulation
    sbatch run_pmf_$LBOUND.sh
    rm run_pmf_$LBOUND.sh
done