# Make a new directory for the simulation
export TAG="poly_sequ_${SEQUNAME}_nchain_${NCHAIN}_temp_${TEMP}_pmf"
mkdir -p "$CWD_PATH/data/$TAG"
cd "$CWD_PATH/data/$TAG"

# Copy the LAMMPS input file to the new directory
cp "$CWD_PATH/parameters/2chain_pmf.in" "$CWD_PATH/data/$TAG/2chain_pmf.in"
cp "$CWD_PATH/parameters/input.colvars" "$CWD_PATH/data/$TAG/input.colvars"
cp "$CWD_PATH/parameters/groups.ndx" "$CWD_PATH/data/$TAG/groups.ndx"

cp "$CWD_PATH/parameters/poly_init.py" "$CWD_PATH/data/$TAG/poly_init.py"
# Initialise box
sed -i "s/SEQUENCE/$sequence/g" poly_init.py
sed -i "s/NCHAIN/$NCHAIN/g" poly_init.py
sed -i "s/LENGTH/$LENGTH/g" poly_init.py

python3 poly_init.py

# Modify the LAMMPS input file with the input parameters
sed -i "s/TEMP/$TEMP/g" 2chain_pmf.in

# Run the LAMMPS simulation
if [ "${NGPU}" == "0" ]
then
    mpirun -np $NCPU $LAMMPS_PATH -in 2chain_pmf.in > poly.out
else
    mpirun -np $NCPU $LAMMPS_PATH -sf gpu -pk gpu $NGPU -in 2chain_pmf.in > poly.out
fi

echo "SIMULATION COMPLETED SUCCESSFULLY"

cd $CWD_PATH
