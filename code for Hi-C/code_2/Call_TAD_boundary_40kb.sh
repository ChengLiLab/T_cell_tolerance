#!/usr/bin/env bash
mkdir TAD_boundary
cd TAD_boundary

for((i=1;i<=23;i++));

do

	matrix2insulation.pl -i ../TAD_files/chr${i}_chr${i}_40k_normalmatrix.txt.gz -is 1000000 -ids 240000 -im mean  -nt 0.1 -v
	  
done

