
###========================================
!/bin/bash
BSUB -n 16
BSUB -R "span[ptile=16]"
BSUB -o lsf.out
BSUB -e lsf.err
BSUB -q "windfall"
BSUB -J GalIC_test
#---------------------------------------------------------------------


module load openmpi
mpirun -np 16 ./../GalIC > myparameterfile.param


###end of script
