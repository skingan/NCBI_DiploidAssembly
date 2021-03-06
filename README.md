# Diploid Assembly Submission to NCBI
NCBI now accepts diploid genome submissions! This repository contains a collection of scripts to aid you in generating the placement file required by NCBI. Details can be found [here](https://www.ncbi.nlm.nih.gov/assembly/docs/submission/). Refer to the [_Arabidopsis thaliana assembly_](https://www.ncbi.nlm.nih.gov/assembly/GCA_001753755.2) from [Chin et al. 2016](https://www.ncbi.nlm.nih.gov/pubmed/27749838) as you prepare your FALCON-Unzip assembly for submission to NCBI.


# Placement File
The placement file for an unscaffolded assembly has the following fields:

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

## Example Placement File
Refer to [AthalPlacementFile.txt](https://github.com/skingan/NCBI_DiploidAssembly/blob/master/AthalPlacementFile.txt) in this repository for the _Arabidpsis thaliana_ file.

# Scripts to generate placement files
Either of the following scripts can be used to generate a placement file.

## nucmer2ncbiPlacement.py
Python script written by [Sarah Kingan](https://github.com/skingan) that generates placement file for diploid genome assembly. This script operates on a directory containing "coords" files created from nucmer alignments between each haplotig and its associated primary contig. The commands to generate "coords" files can be found within script on lines 23-43.

## generate_placement.py
python script written by [Jason Chin](https://github.com/pb-jchin) that performs nucmer alignments, generates coords file, and created placement file for NCBI submission. Dependencies include [FastaReader.py](https://github.com/PacificBiosciences/FALCON/blob/master/falcon_kit/FastaReader.py) and [mummer](http://mummer.sourceforge.net/).



THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
