# Make a new directory for the simulation
export TAG="poly_sequ_${SEQUNAME}_nchain_${NCHAIN}_znet_${ZNET}_nsalt_${NSALT}_temp_${TEMP}_pmf"
mkdir -p "$CWD_PATH/data/$TAG/$LBOUND"
cd "${CWD_PATH}/data/${TAG}/${LBOUND}"


if [ -f poly.out ]
then 
    if grep -Fq "Total wall time:" poly.out
    then
        echo "SIMULATION COMPLETED SUCCESSFULLY"
    else
        echo "SIMULATION STARTED BUT DIDNT COMPLETE"
        cp "$CWD_PATH/parameters/2chain_pmf_restart.in" "$CWD_PATH/data/$TAG/$LBOUND/2chain_pmf_restart.in"

        sed -i "s/TEMP/$TEMP/g" 2chain_pmf_restart.in
        sed -i "s/STEPS/$STEPS/g" 2chain_pmf_restart.in

        if test -f "restart2"
        then
            if [ "restart1" -nt "restart2" ]
            then
                sed -i "s/restartN/restart1/g" 2chain_pmf_restart.in
            else
                sed -i "s/restartN/restart2/g" 2chain_pmf_restart.in
            fi
        else
            sed -i "s/restartN/restart1/g" 2chain_pmf_restart.in
        fi



        if [ "${NGPU}" == "0" ]
        then
            mpirun -np $NCPU $LAMMPS_PATH -in 2chain_pmf_restart.in >> poly.out
        else
            mpirun -np $NCPU $LAMMPS_PATH -sf gpu -pk gpu $NGPU -in 2chain_pmf_restart.in >> poly.out
        fi

    fi
else
    # Copy the LAMMPS input file to the new directory
    cp "$CWD_PATH/parameters/2chain_pmf.in" "$CWD_PATH/data/$TAG/$LBOUND/2chain_pmf.in"
    cp "$CWD_PATH/parameters/input.colvars" "$CWD_PATH/data/$TAG/$LBOUND/input.colvars"
    cp "$CWD_PATH/parameters/groups.ndx" "$CWD_PATH/data/$TAG/$LBOUND/groups.ndx"

    cp "$CWD_PATH/parameters/poly_init_pmf.py" "$CWD_PATH/data/$TAG/$LBOUND/poly_init.py"
    # Initialise box
    sed -i "s/SEQUENCE/$sequence/g" poly_init.py
    sed -i "s/NCHAIN/$NCHAIN/g" poly_init.py
    sed -i "s/NSALT/$NSALT/g" poly_init.py
    sed -i "s/LENGTHX/$LENGTHX/g" poly_init.py
    sed -i "s/LENGTHY/$LENGTHY/g" poly_init.py
    sed -i "s/LENGTHZ/$LENGTHZ/g" poly_init.py
    sed -i "s/LBOUND/$LBOUND/g" poly_init.py
    sed -i "s/UBOUND/$UBOUND/g" poly_init.py

    python3 poly_init.py

    # Modify the LAMMPS input file with the input parameters
    sed -i "s/TEMP/$TEMP/g" 2chain_pmf.in
    sed -i "s/STEPS/$STEPS/g" 2chain_pmf.in

    sed -i "s/LBOUND/$LBOUND/g" input.colvars
    sed -i "s/UBOUND/$UBOUND/g" input.colvars

    # Run the LAMMPS simulation
    if [ "${NGPU}" == "0" ]
    then
        mpirun -np $NCPU $LAMMPS_PATH -in 2chain_pmf.in > poly.out
    else
        mpirun -np $NCPU $LAMMPS_PATH -sf gpu -pk gpu $NGPU -in 2chain_pmf.in > poly.out
    fi

    echo "SIMULATION COMPLETED SUCCESSFULLY"
fi

cd $CWD_PATH
