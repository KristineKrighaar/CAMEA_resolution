#!/bin/bash

    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_old_0_v4_10 -n 1000000000 slitwidth=0.03 SourceE=3.2  DeltaSourceE=0.3  EI=3.2  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0 
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_old_1_v4_10 -n 1000000000 slitwidth=0.03 SourceE=3.38  DeltaSourceE=0.3  EI=3.38  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0   
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_old_2_v4_10 -n 1000000000 slitwidth=0.03 SourceE=3.58  DeltaSourceE=0.3  EI=3.58  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_old_3_v4_10 -n 1000000000 slitwidth=0.03 SourceE=3.79  DeltaSourceE=0.3  EI=3.79  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_old_4_v4_10 -n 1000000000 slitwidth=0.03 SourceE=4.04  DeltaSourceE=0.3  EI=4.04  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_old_5_v4_10 -n 1000000000 slitwidth=0.03 SourceE=4.32  DeltaSourceE=0.3  EI=4.32  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_old_6_v4_10 -n 1000000000 slitwidth=0.03 SourceE=4.64  DeltaSourceE=0.3  EI=4.64  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_old_7_v4_10 -n 1000000000 slitwidth=0.03 SourceE=5  DeltaSourceE=0.3  EI=5  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0  Be_filter=1  old_mono=1  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0
    
    

    
    
    
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_new_0_v4_10 -n 1000000000 slitwidth=0.03 SourceE=3.2  DeltaSourceE=0.3  EI=3.2  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0.5  Be_filter=1  old_mono=0  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0 
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_new_1_v4_10 -n 1000000000 slitwidth=0.03 SourceE=3.38  DeltaSourceE=0.3  EI=3.38  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0.5  Be_filter=1  old_mono=0  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0   
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_new_2_v4_10 -n 1000000000 slitwidth=0.03 SourceE=3.58  DeltaSourceE=0.3  EI=3.58  A3=0  A4=-90  SAMPLE=3 file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0.5  Be_filter=1  old_mono=0  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_new_3_v4_10 -n 1000000000 slitwidth=0.03 SourceE=3.79  DeltaSourceE=0.3  EI=3.79  A3=0  A4=-90  SAMPLE=3 file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0.5  Be_filter=1  old_mono=0  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_new_4_v4_10 -n 1000000000 slitwidth=0.03 SourceE=4.04  DeltaSourceE=0.3  EI=4.04  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0.5  Be_filter=1  old_mono=0  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_new_5_v4_10 -n 1000000000 slitwidth=0.03 SourceE=4.32  DeltaSourceE=0.3  EI=4.32  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0.5  Be_filter=1  old_mono=0  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_new_6_v4_10 -n 1000000000 slitwidth=0.03 SourceE=4.64  DeltaSourceE=0.3  EI=4.64  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0.5  Be_filter=1  old_mono=0  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0
    
    mcrun -c --mpi=8 FullInstrument_v4_elastic_resolution.instr -d res_scan_new_7_v4_10 -n 1000000000 slitwidth=0.03 SourceE=5  DeltaSourceE=0.3  EI=5  A3=0  A4=-90  SAMPLE=3  file_name=EGECE_13WaveLmax1_geometry.dat.txt  RV_mono_Bool=0.5  RH_mono_Bool=0.5  Be_filter=1  old_mono=0  sampleHeight=0.0105  sampleRadius=0.005  samplePosY=0

