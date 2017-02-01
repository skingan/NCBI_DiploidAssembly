# Diploid Assembly Submission to NCBI
NCBI now accepts diploid genome submissions! Refer to the [Arabidopsis thaliana assembly] (https://www.ncbi.nlm.nih.gov/assembly/GCA_001753755.2) from [Chin et al. 2016] (https://www.ncbi.nlm.nih.gov/pubmed/27749838) as you prepare your FALCON-Unzip assembly for submission to NCBI.

# Placement File
This repositor contains a collection of scripts to aid you in generating the placement file required by NCBI. The placement file for an unscaffolded assembly has the following fields:

1. alt_asm_name: name of the assembly-unit that includes the alternate scaffold.
2. prim_asm_name: name of the assembly-unit on which the alternate scaffold is being placed. Expected to be 'Primary Assembly' in most cases.
3. alt_scaf_name: name of the alternate scaffold being placed
4. parent_type: type of object on which the alternate scaffold is being placed, either CHROMOSOME or SCAFFOLD
5. parent_name: name of the object on which the alternate scaffold is being placed (can be either a chromosome or a scaffold)
6. ori: orientation of the alignment, '+', '-', 'b'(mixed)
7. alt_scaf_start: start of the placement on the alternate scaffold (in 1 base coordinates)
8. alt_scaf_stop: end of the placement on the alternate scaffold (in 1 base coordinates)
9. parent_start: start of the placement on the parent sequence (in 1 base coordinates)
10. parent_stop: end of the placement on the parent sequence (in 1 base coordinates)
11. alt_start_tail: number of bases at the start of the alternate scaffold not involved in the placement
12. alt_stop_tail: number of bases at the end of the alternate scaffold not involved in the placement

# Example Placement File
Refer to AthalPlacementFile.txt in this repository.

# nucmer2ncbiPlacement.py
python script that generates placement file for diploid genome assembly. This script operates on a directory containing "coords" files created from nucmer alignments between each haplotig and its associated primary contig. The commands to generate "coords" files can be found within script as commented lines.

# 
