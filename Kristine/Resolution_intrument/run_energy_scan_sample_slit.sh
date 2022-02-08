#!/bin/bash

mcrun FullInstrument_v4.7_frontend.instr -d vir_slit_optimal_scan_1 -n 100000000 -N 11 slitwidth=0.01,0.1  virtual_slit_d=0.38  SourceE=5  DeltaSourceE=0.3  EI=5  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0.5  Be_filter=1  old_mono=0  sampleHeight=0.01  sampleRadius=0.005  samplePosY=0




    mcrun -c --mpi=4 FullInstrument_v4_elastic_resolution.instr -d /home/kristine/Dropbox/Masters_2/CAMEA_Resolution_project/Resolution_intrument/energy_scan_old_0 -n 100000000 SourceE=3.2  DeltaSourceE=1  EI=3.2  A3=0  A4=-90  SAMPLE=4  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.01  sampleRadius=0.005  samplePosY=0 
    
    mcrun -c --mpi=4 FullInstrument_v4_elastic_resolution.instr -d /home/kristine/Dropbox/Masters_2/CAMEA_Resolution_project/Resolution_intrument/energy_scan_old_1 -n 100000000 SourceE=3.38  DeltaSourceE=1  EI=3.38  A3=0  A4=-90  SAMPLE=4  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.01  sampleRadius=0.005  samplePosY=0   
    
    mcrun -c --mpi=4 FullInstrument_v4_elastic_resolution.instr -d /home/kristine/Dropbox/Masters_2/CAMEA_Resolution_project/Resolution_intrument/energy_scan_old_2 -n 100000000 SourceE=3.58  DeltaSourceE=1  EI=3.58  A3=0  A4=-90  SAMPLE=4  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.01  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=4 FullInstrument_v4_elastic_resolution.instr -d /home/kristine/Dropbox/Masters_2/CAMEA_Resolution_project/Resolution_intrument/energy_scan_old_3 -n 100000000 SourceE=3.79  DeltaSourceE=1  EI=3.79  A3=0  A4=-90  SAMPLE=4  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.01  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=4 FullInstrument_v4_elastic_resolution.instr -d /home/kristine/Dropbox/Masters_2/CAMEA_Resolution_project/Resolution_intrument/energy_scan_old_4 -n 100000000 SourceE=4.04  DeltaSourceE=1  EI=4.04  A3=0  A4=-90  SAMPLE=4  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.01  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=4 FullInstrument_v4_elastic_resolution.instr -d /home/kristine/Dropbox/Masters_2/CAMEA_Resolution_project/Resolution_intrument/energy_scan_old_5 -n 100000000 SourceE=4.32  DeltaSourceE=1  EI=4.32  A3=0  A4=-90  SAMPLE=4  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.01  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=4 FullInstrument_v4_elastic_resolution.instr -d /home/kristine/Dropbox/Masters_2/CAMEA_Resolution_project/Resolution_intrument/energy_scan_old_6 -n 100000000 SourceE=4.64  DeltaSourceE=1  EI=4.64  A3=0  A4=-90  SAMPLE=4  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.01  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=4 FullInstrument_v4_elastic_resolution.instr -d /home/kristine/Dropbox/Masters_2/CAMEA_Resolution_project/Resolution_intrument/energy_scan_old_7 -n 100000000 SourceE=5  DeltaSourceE=1  EI=5  A3=0  A4=-90  SAMPLE=4  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.01  sampleRadius=0.005  samplePosY=0
    
    
    
   

