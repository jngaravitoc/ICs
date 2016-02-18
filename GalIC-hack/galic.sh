###========================================
#!/bin/bash
#BSUB -n 32
#BSUB -R "span[ptile=16]"
#BSUB -o lsf.out
#BSUB -e lsf.err
#BSUB -q "medium"
#BSUB -J GalIC_test
#---------------------------------------------------------------------

module load gsl
#module load openmpi
mpirun -np 32 ./GalIC  MWbesla07.param


###end of script
