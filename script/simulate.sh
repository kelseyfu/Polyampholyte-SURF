# Make a new directory for the simulation
export TAG="poly_sequ_${SEQUNAME}_nchain_${NCHAIN}_temp_${TEMP}"
mkdir -p "$CWD_PATH/data/$TAG"

# Copy the LAMMPS input file to the new directory
cp "$CWD_PATH/parameters/poly.in" "$CWD_PATH/data/$TAG/poly.in"
cp "$CWD_PATH/parameters/poly_init.py" "$CWD_PATH/data/$TAG/poly_init.py"
cd "$CWD_PATH/data/$TAG"

# Initialise box
sed -i '' "s/SEQUENCE/$sequence/g" poly_init.py
sed -i '' "s/NCHAIN/$NCHAIN/g" poly_init.py

python3 poly_init.py

# Modify the LAMMPS input file with the input parameters
sed -i '' "s/TEMP/$TEMP/g" poly.in

# Run the LAMMPS simulation
$LAMMPS_PATH -in poly.in > poly.out

cd $CWD_PATH
