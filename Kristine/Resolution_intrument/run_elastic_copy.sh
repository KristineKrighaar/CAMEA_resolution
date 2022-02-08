#!/bin/bash
    
mcrun -c --mpi=4 FullInstrument_v4.5_source_sample_pos.instr -d test_energi_monitors_5 -n 1000000000 SourceE=4.2  DeltaSourceE=2  EI=4.2  A3=0  A4=-90  file_name=EGECE_13WaveLmax1_geometry.dat.txt  Be_filter=1  sampleHeight=0.01  sampleRadius=0.005  samplePosY=0

