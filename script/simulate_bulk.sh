# Make a new directory for the simulation
export TAG="poly_sequ_${SEQUNAME}_nchain_${NCHAIN}_znet_${ZNET}_nsalt_${NSALT}_temp_${TEMP}_bulk"
mkdir -p "$CWD_PATH/data/$TAG"
cd "$CWD_PATH/data/$TAG"

if [ -f poly.out ]
then 
    if grep -Fq "Total wall time:" poly.out
    then
        echo "SIMULATION COMPLETED SUCCESSFULLY"
    else
        echo "SIMULATION STARTED BUT DIDNT COMPLETE"
        cp "$CWD_PATH/parameters/poly_bulk_restart.in" "$CWD_PATH/data/$TAG/poly_bulk_restart.in"

        sed -i "s/TEMP/$TEMP/g" poly_bulk_restart.in
        sed -i "s/STEPS/$STEPS/g" poly_bulk_restart.in

        if [ "${NGPU}" == "0" ]
        then
            mpirun -np $NCPU $LAMMPS_PATH -in poly_bulk_restart.in >> poly.out
        else
            mpirun -np $NCPU $LAMMPS_PATH -sf gpu -pk gpu $NGPU -in poly_bulk_restart.in >> poly.out
        fi

    fi
else
    # Copy the LAMMPS input file to the new directory
    cp "$CWD_PATH/parameters/poly_bulk.in" "$CWD_PATH/data/$TAG/poly_bulk.in"
    cp "$CWD_PATH/parameters/poly_init.py" "$CWD_PATH/data/$TAG/poly_init.py"
    # Initialise box
    sed -i "s/SEQUENCE/$sequence/g" poly_init.py
    sed -i "s/NCHAIN/$NCHAIN/g" poly_init.py
    sed -i "s/NSALT/$NSALT/g" poly_init.py
    sed -i "s/LENGTHX/$LENGTHX/g" poly_init.py
    sed -i "s/LENGTHY/$LENGTHY/g" poly_init.py
    sed -i "s/LENGTHZ/$LENGTHZ/g" poly_init.py

    python3 poly_init.py

    # Modify the LAMMPS input file with the input parameters
    sed -i "s/TEMP/$TEMP/g" poly_bulk.in
    sed -i "s/STEPS/$STEPS/g" poly_bulk.in

    # Run the LAMMPS simulation
    if [ "${NGPU}" == "0" ]
    then
        mpirun -np $NCPU $LAMMPS_PATH -in poly_bulk.in > poly.out
    else
        mpirun -np $NCPU $LAMMPS_PATH -sf gpu -pk gpu $NGPU -in poly_bulk.in > poly.out
    fi

    echo "SIMULATION COMPLETED SUCCESSFULLY"
fi

cd $CWD_PATH
