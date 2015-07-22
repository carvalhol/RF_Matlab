#!/bin/bash

clear


for i in {1..1}
do
   echo "Running $i"
cd 3D_W
cd ISO-i_3D_W
cd 1p_L-1_ISO-i_3D_W
qsub run1.pbs
cd ../../ISO-i_3D_W
cd 2p_L-2_ISO-i_3D_W
qsub run2.pbs
cd ../../ISO-i_3D_W
cd 4p_L-3_ISO-i_3D_W
qsub run4.pbs
cd ../../ISO-i_3D_W
cd 8p_L-4_ISO-i_3D_W
qsub run8.pbs
cd ../../ISO-i_3D_W
cd 16p_L-5_ISO-i_3D_W
qsub run16.pbs
cd ../../ISO-i_3D_W
cd 32p_L-6_ISO-i_3D_W
qsub run32.pbs
cd ../../ISO-i_3D_W
cd 64p_L-7_ISO-i_3D_W
qsub run64.pbs
cd ../../ISO-i_3D_W
cd 128p_L-8_ISO-i_3D_W
qsub run128.pbs
cd ../../ISO-i_3D_W
cd 256p_L-9_ISO-i_3D_W
qsub run256.pbs
cd ../../ISO-i_3D_W
cd 512p_L-10_ISO-i_3D_W
qsub run512.pbs
cd ../../SHINO_3D_W
cd 1p_L-1_SHINO_3D_W
qsub run1.pbs
cd ../../SHINO_3D_W
cd 2p_L-2_SHINO_3D_W
qsub run2.pbs
cd ../../SHINO_3D_W
cd 4p_L-3_SHINO_3D_W
qsub run4.pbs
cd ../../SHINO_3D_W
cd 8p_L-4_SHINO_3D_W
qsub run8.pbs
cd ../../SHINO_3D_W
cd 16p_L-5_SHINO_3D_W
qsub run16.pbs
cd ../../SHINO_3D_W
cd 32p_L-6_SHINO_3D_W
qsub run32.pbs
cd ../../SHINO_3D_W
cd 64p_L-7_SHINO_3D_W
qsub run64.pbs
cd ../../SHINO_3D_W
cd 128p_L-8_SHINO_3D_W
qsub run128.pbs
cd ../../SHINO_3D_W
cd 256p_L-9_SHINO_3D_W
qsub run256.pbs
cd ../../SHINO_3D_W
cd 512p_L-10_SHINO_3D_W
qsub run512.pbs
cd ../../RANDO_3D_W
cd 1p_L-1_RANDO_3D_W
qsub run1.pbs
cd ../../RANDO_3D_W
cd 2p_L-2_RANDO_3D_W
qsub run2.pbs
cd ../../RANDO_3D_W
cd 4p_L-3_RANDO_3D_W
qsub run4.pbs
cd ../../RANDO_3D_W
cd 8p_L-4_RANDO_3D_W
qsub run8.pbs
cd ../../RANDO_3D_W
cd 16p_L-5_RANDO_3D_W
qsub run16.pbs
cd ../../RANDO_3D_W
cd 32p_L-6_RANDO_3D_W
qsub run32.pbs
cd ../../RANDO_3D_W
cd 64p_L-7_RANDO_3D_W
qsub run64.pbs
cd ../../RANDO_3D_W
cd 128p_L-8_RANDO_3D_W
qsub run128.pbs
cd ../../RANDO_3D_W
cd 256p_L-9_RANDO_3D_W
qsub run256.pbs
cd ../../RANDO_3D_W
cd 512p_L-10_RANDO_3D_W
qsub run512.pbs
cd ../../../
sleep 1
done

qstat -u carvalhol
