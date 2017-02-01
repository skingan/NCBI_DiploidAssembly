# Jason Chin
# Pacific Biosciences
# 2016

from falcon_kit.FastaReader import FastaReader
import os


primary_contigs = {}
f = FastaReader("cns_p_ctg.fasta")
for r in f:
    rname = r.name.split("|")[0]
    primary_contigs.setdefault( rname, (r.sequence, []))
primary_contigs.setdefault( "NA", ("", []))

f = FastaReader("cns_h_ctg.fasta")
all_h_ctg = set()
for r in f:
    rname = r.name.split("|")[0]
    p_name = rname.split("_")[0]
    all_h_ctg.add(rname)
    if p_name in primary_contigs:
        primary_contigs[p_name][1].append( (rname, r.sequence) )
        #print rname, p_name
    else:
        primary_contigs["NA"][1].append( (rname, r.sequence) )


data = []
place_h_ctg = set()
for ctg in primary_contigs:
    if ctg == "NA":
        continue
    p_ctg_seq = primary_contigs[ctg][0]
    if len(primary_contigs[ctg][1]) == 0:
        continue
    
    with open("%s.fa" % ctg, "w") as f:
        print >> f, ">%s" % ctg
        print >> f, p_ctg_seq

    with open("h_%s.fa" % ctg, "w") as f:
        for h_ctg, h_ctg_seq in primary_contigs[ctg][1]:
            print >>f, ">%s" % h_ctg
            print >>f,  h_ctg_seq

    os.system("~/build/MUMmer3.23/nucmer -mum %s.fa h_%s.fa -p %s" % (ctg, ctg, ctg))
    os.system("~/build/MUMmer3.23/delta-filter -g -q  %s.delta > %s_1.delta" % (ctg, ctg))

    os.system("~/build/MUMmer3.23/show-coords -r -H -L 500 -l -T %s_1.delta > %s_1.coor" % (ctg, ctg))

    hctg_data = {}
    with open("%s_1.coor" % ctg) as f:
        for row in f:
            row = row.strip().split()
            s1, e1, s2, e2, l1, l2, idt, lr, lq, p_ctg, h_ctg = row
            s1 = int(s1)
            e1 = int(e1)
            s2 = int(s2)
            e2 = int(e2)
            l1 = int(l1)
            l2 = int(l2)
            lr = int(lr)
            lq = int(lq)
            idt = float(idt)
            if idt < 85:
                continue
            #if l1 < 2000:
            #    continue
            hctg_data.setdefault(h_ctg, [])
            hctg_data[ h_ctg ].append( (p_ctg, s1, e1, s2, e2, l1, l2, lr, lq, idt) ) 

    for h_ctg in hctg_data:
        d = hctg_data[ h_ctg ]
        d.sort()
        p_ctg, s1, e1, s2, e2, l1, l2, lr, lq, idt = d[0]
        parent_start = s1
        start_tail = s2 - 1
        scaf_start = s2
        p_ctg, s1, e1, s2, e2, l1, l2, lr, lq, idt = d[-1]
        parent_stop = e1
        stop_tail = lq - e2 
        scaf_stop = e2

        ori = "+"
        if scaf_start > scaf_stop:
            scaf_start, scaf_stop = scaf_stop, scaf_start
            start_tail, stop_tail = lq - stop_tail - 1, lq - start_tail - 1
            ori = "-"
        place_h_ctg.add(h_ctg)
        data.append( (h_ctg, "SCAFFOLD", p_ctg, ori, scaf_start, scaf_stop, parent_start, parent_stop, start_tail, stop_tail) )

"""
#asm_id_from    asm_id_to       scaf_name_from  obj_type_to     obj_name_to     ori     scaf_start      scaf_stop       parent_start parent_stop start_tail stop_tail
haplotigs       Primary Assembly                haplotig1       SCAFFOLD        contig1 +               1               439505       88891175    89330679   0 0
haplotigs
"""
data.sort(key=lambda x: (x[2],x[6]))
print """#asm_id_from    asm_id_to       scaf_name_from  obj_type_to     obj_name_to     ori     scaf_start      scaf_stop       parent_start parent_stop start_tail stop_tail"""
for d in data:
    s= ["haplotigs", "Primary Assembly"] + [str(c) for c in d]
    print "\t".join(s)
