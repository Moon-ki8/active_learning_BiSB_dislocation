#!/bin/bash
#SBATCH --job-name=main_train
#SBATCH --output=result.out
#SBATCH --error=error.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --partition=cpu
#SBATCH --account=bcol-delta-cpu
#SBATCH --time=10:00:00
#SBATCH --mail-user=choi0652@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=0  # This requests all available memory on the node

module load gcc/11.4.0
module load openmpi/4.1.6
module load openblas/0.3.25
module load netlib-lapack
module load intel-oneapi-mkl/2024.0.0

export QUIP_ARCH=linux_x86_64_gfortran
export LD_LIBRARY_PATH=/sw/spack/deltas11-2023-03/apps/linux-rhel8-zen3/gcc-11.4.0/netlib-lapack-3.11.0-dgkcd33/lib64:$LD_LIBRARY_PATH
source /u/mchoi3/myenv/bin/activate
export OMP_NUM_THREADS=128

module load python/3

## INPUT PARAMETER 
pot_index=4 # 0 if it is initial
num_proc=256
export cur_path=$(pwd)
export pot_dir='/scratch/bdai/mchoi3/Active_learning_BiSb_v7_force_error/POTENTIAL'

rm fnorm.txt
rm error_max.txt

#### FUNCTIONS ####
# Evaluation of current data file (after.data)
EVAL(){
cd ${cur_path}
cp after_$2.data after.data
python utils.py "data2xyz"
/u/mchoi3/QUIP2/QUIP/build/linux_x86_64_gfortran_openmpi+openmp/quip output_file=quip.out verbosity=ANALYSIS atoms_filename=after.xyz calc_args="local_gap_variance" param_filename=${pot_dir}/gap_BiSb_$1.xml E=T F=T L=T > eval.out
python utils.py "eval_gap2xyz"
mv eval.out EVAL/eval_$1.out
cp eval.xyz EVAL/eval.xyz
cp eval.xyz EVAL/eval_$1.xyz
max_value=$(python utils.py "find_max_gap_value") # Maximum error value 
}

# Make a representative image of atoms near the maximum error
MAKEREP(){
pot_index=$((pot_index+1))
python utils.py "create_rep_new"
mv rep.xyz DFT/rep.xyz
}

MIN(){
cd ${pot_dir}
var="gap_BiSb_$pot_index.xml.sparseX.GAP"
FILES=$(find . -type f -name "*$var*")
for FILE in $FILES; do
  extracted=$(echo "$FILE" | sed 's/.*\.\(GAP_[0-9_\-]*\)[0-9]/\1/') # Potential file name extraction
done
cd ${cur_path}

# Write LAMMPS min input script 
echo "units metal" > min.in
echo "boundary f p f " >> min.in
echo "read_data after_$1.data" >> min.in
echo "pair_style quip" >> min.in
echo "pair_coeff * * ${pot_dir}/gap_BiSb_${pot_index}.xml \"Potential xml_label=${extracted}\" 83 51" >> min.in
echo "thermo 100" >> min.in
echo "thermo_style    custom step etotal fmax fnorm pxx pzz pxz" >> min.in 
echo "thermo_modify   format line \"%d %18.17e %18.17e %18.17e %18.17e %18.17e %18.17e\"" >> min.in
echo "mass 1 208.9804" >> min.in
echo "mass 2 121.76" >> min.in
echo "dump           dump1 all custom 100 displaced.dump id type x y z fx fy fz" >> min.in
echo "dump_modify    dump1 format line \"%d %d %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\"" >> min.in
# Constrain 
echo "variable L equal lx/2" >> min.in
echo "variable rin equal \$L-8" >> min.in
echo "variable rout equal \$L" >> min.in
echo "region cout cylinder y \$L \$L \${rout} EDGE EDGE" >> min.in
echo "region cin cylinder y \$L \$L \${rin} EDGE EDGE" >> min.in
echo "group outatom region cout" >> min.in
echo "group inatom region cin" >> min.in
echo "group cedge subtract outatom inatom" >> min.in
echo "fix 2 cedge setforce 0.0 NULL 0.0" >> min.in 
echo "min_style      fire" >> min.in
echo "min_modify      integrator verlet tmax 0.2 tmin 0.0" >> min.in
echo "minimize       0.0 1e-8  300 700000" >> min.in
echo "variable fn equal fnorm" >> min.in
echo "variable minflag equal 0" >> min.in
echo "if \"\${fn} < 1e-8\" then \"variable minflag equal 1\"" >> min.in
echo "print \"convergence=\${minflag}\" file conv_min.txt" >> min.in
echo "write_data after_$2.data" >> min.in
echo "print \"\${fn}\" append fnorm.txt" >> min.in
srun --mpi=pmi2 /u/mchoi3/LAMMPS/lammps-2Aug2023/build/lmp_mpi -in  min.in > out.txt 
file="conv_min.txt"
## Read the converge value into a variable
converge=$(grep 'convergence=' "$file" | cut -d'=' -f2)
}
DFT(){
	cd DFT
    mv out.txt out_prev.txt
	python scf_generator.py
	srun --mpi=pmi2 /u/mchoi3/qe_CPU/qe-7.2/bin/pw.x -in scf.in > out.txt
	python qe2xyz.py
    python config_type_set.py
	cp rep_force.xyz ${pot_dir}/rep_$1.xyz
    cd ${cur_path}	
}

# Update potential using new dft data
UPDATEPOT(){
cd ${pot_dir}
cat rep_$1.xyz >> DB_total.xyz
time /u/mchoi3/QUIP2/QUIP/build/linux_x86_64_gfortran_openmpi+openmp/gap_fit at_file=DB_total.xyz \
gap={ distance_Nb order=2 cutoff=5.0 n_sparse=50 covariance_type=ard_se \
     delta=0.5 theta_uniform=2.5 sparse_method=uniform compact_clusters=True \
     :soap l_max=8 \
                       n_max=8 \
                       cutoff=5.0 \
                       atom_sigma=0.12 \
                       radial_scaling=-0.5 \
                       cutoff_transition_width=0.5 \
                       central_weight=1.0 \
                       config_type_n_sparse = {binary:600:surface:200:rep:500} \
                       delta=0.5 \
                       covariance_type=dot_product \
                       zeta=4 \
                       sparse_method=random \
                  :soap l_max=8 \
                       n_max=8 \
                       cutoff=5.0 \
                       atom_sigma=0.10 \
                       radial_scaling=-0.5 \
                       cutoff_transition_width=0.5 \
                       central_weight=1.0 \
                       delta=0.5 \
                       covariance_type=dot_product \
                       config_type_n_sparse = {binary:600:surface:200:rep:500} \
                       zeta=4 \
                       sparse_method=random} \
default_sigma={0.001 0.001 1.0 0.0} \
config_type_kernel_regularisation={binary:0.001:0.01:1.0:0.0:surface:0.001:0.001:1.0:0.0:rep:0.001:0.001:1.0:0.0} \
energy_parameter_name=energy \
force_parameter_name=forces \
virial_parameter_name=stress \
sparse_jitter=1.0e-8 \
do_copy_at_file=F \
sparse_separate_file=T \
gp_file=gap_BiSb_$1.xml \
openmp_chunk_size=100 > pot_out_$1.txt
cd ${cur_path}
}

########### MAIN ############
runstep=300
threshold=0.015 # Force variacne threshold
init_run=67800
runs=$init_run
count=0
cd ${cur_path}

# Loop over a range of numbers
for num in {1..4}; do
    echo "num= $num"
    flag_high_error=0
    minflag=false
    # Inner loop to handle evaluation and minimization steps
    while [ "$flag_high_error" -eq 0 ]; do
        # Perform minimization if count is greater than 0
        if [ "$minflag" = true ]; then
            prevruns=$((runstep * (count - 1) + init_run))
            runs=$((runstep * count + init_run))
            echo "Run Minimization: $runs steps"
            MIN $prevruns $runs
	        minflag=false
            mv log.lammps log_$runs.lammps
        fi
        # Evaluate potential energy
        EVAL $pot_index $runs

        # Check maximum uncertainty error
        echo "Max uncertainty error: $max_value"
        echo "$runs $max_value" >> error_max.txt # Record error_max
        compare=$(echo "$max_value > $threshold" | bc)
        if [ "$compare" -eq 1 ]; then
            echo "Max error exceeds threshold, create rep image"
            MAKEREP
            echo "pot_index = ${pot_index}"
            DFT $pot_index
            flag_high_error=1
            UPDATEPOT $pot_index
            echo "Update done"
	    else
            echo "Max error does not exceed threshold, continuing minimization"
            minflag=true
            # Increment count for the next step
            count=$((count + 1))
        fi

    done
done

cp result.out result_4.txt

