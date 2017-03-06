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
import numpy as np
import subprocess

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

def file_len(fname):
	#http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
	p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	result, err = p.communicate()
	if p.returncode != 0:
		raise IOError(err)
	return int(result.strip().split()[0])



fileList = glob.glob('*.coords')
for file in fileList:
# htigs that align to primary
	if file_len(file) >= 5:
		output=['haplotigs','Primary_Assembly', '000000F_000', 'SCAFFOLD', '000000F', '+', 'hstart', 'hstop', 'pstart', 'pstop', '0', '0']
		d=np.loadtxt(open(file, "rb"), skiprows=4, dtype="str", ndmin=2)
		output[2]=d[0,10]	#htig name
		output[4]=d[0,9]	#pcontig name
		h=np.array(d[:,[2,3]],dtype=int) # coax into right data type to use amin
		p=np.array(d[:,[0,1]],dtype=int)
		output[6]=np.amin(h) # hstart
		output[7]=np.amax(h) # hstop
		output[8]=np.amin(p) # pstart
		output[9]=np.amax(p) # pstop
# tail regions
		if int(output[6]) != 1: 		# if aln doesn't start at begining of htig
			output[10]=int(output[6])-1
		if int(output[7]) != int(d[0,8]): 	# if aln doesn't stop at end of htig
			output[11]=int(d[0,8])-int(output[7])
# alignment orientation
		if np.argmin(p) > np.argmax(p):
			output[5] = '-'
# htigs that don't align to primary
	else:
		output=['haplotigs','Primary Assembly',file.split('.')[0],'na','na','na','na','na','na','na','na','na']
# print to stdout
	print '\t'.join(str(o) for o in output)

