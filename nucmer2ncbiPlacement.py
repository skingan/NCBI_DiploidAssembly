#!/usr/bin/env python

# Sarah B. Kingan
# 28 January 2017
# 
# Pacific Biosciences
# Applications Lab
#
# Convert directory of nucmer.coords files for haplotigs aligned
# to primary contiga into ncbi placement file
#
###########################################################################################################

# import libraries
import glob

###########################################################################################################

# print header
header=['#alt_asm_name','prim_asm_name','alt_scaf_name','parent_type','parent_name','ori','alt_scaf_start','alt_scaf_stop','parent_start','parent_stop','alt_start_tail','alt_stop_tail']
print '\t'.join(str(h) for h in header)

# process each coords file
# coords files created with shell script:
#
#REF_DIR=/path/to/dir/containing/reference/file
#REF=$REF_DIR/referenceFilePrimaryContigsAndHaplotigs.fasta
#HAP_LIST=`cat $REF_DIR/listOfHaplotigIDs.txt`
#
#
#module load samtools
#module load mummer
#
#for h in $HAP_LIST
#do
#	samtools faidx $REF $h | sed 's/|quiver|arrow|arrow//g' > hap.fa  	# query file of haplotig seq
#	p=`echo $h | sed 's/_[0-9]\{3\}//'`					# associated primary contig id
#	samtools faidx $REF $p | sed 's/|quiver|arrow|arrow//g' > prim.fa	# reference file of primary contig seq
#	prefix=`echo $h | sed 's/|quiver|arrow|arrow//'`			# prefix for nucmer run
#	nucmer -maxmatch prim.fa hap.fa -prefix $prefix				# nucmer run
#	delta-filter -g $prefix.delta > $prefix.global.delta			# 1-1 global aln without rearragments
#	show-coords -qTl $prefix.global.delta > $prefix.coords			# coords format, tab delimited, sorted by query (htig), with ref and query lengths
#done;

fileList = glob.glob('*.coords')
for file in fileList:
	output=['haplotigs','Primary Assembly', '000000F_000', 'SCAFFOLD', '000000F', '+', 'hstart', 'hstop', 'pstart', 'pstop', '0', '0']
	with open(file, 'r') as fin:
		data=[line.split() for line in fin]
# haplotigs that align to primary contig
		if len(data) > 4:
			output[2]=data[4][10]	#htig name
			output[4]=data[4][9]	#pcontig name
			output[6]=data[4][2]	#hstart
			output[7]=data[-1][3]	#hstop
			output[8]=data[4][0]	#pstart
			output[9]=data[-1][1]	#pstop
# tail regions
			if int(output[6]) != 1: # if aln doesn't start at 1
				output[10]=int(output[6])-1 #hstart tail

			if int(output[7]) != int(data[4][8]): # if aln doesn't end at length
				output[11]=int(data[4][8])-int(output[7]) #hstop tail
# haplotigs without alignment to primary
		else:
			output=['haplotigs','Primary Assembly',file.split('.')[0],'na','na','na','na','na','na','na','na','na']
# print output for each haplotig
		print '\t'.join(str(o) for o in output)
