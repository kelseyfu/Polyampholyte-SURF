# Make a new directory for the simulation
export TAG="poly_sequ_${SEQUNAME}_nchain_${NCHAIN}_temp_${TEMP}"
echo $TAG
cd "$CWD_PATH/data/$TAG"
echo $(pwd)

# Copy the LAMMPS analysis file to the new directory
cp "$CWD_PATH/parameters/poly_analysis.in" "$CWD_PATH/data/$TAG/poly_analysis.in"
echo $(ls)

# Run the LAMMPS simulation
$LAMMPS_PATH -in poly_analysis.in > poly_analysis.out

cd $CWD_PATH
