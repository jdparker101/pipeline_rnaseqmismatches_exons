'''
count_mismatches.py - count the number of mismatches per gene
====================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Count the number of high quality mismatches per gene per base sequenced. Will discard reads marked as duplicate

Usage
-----

.. Example use case

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
from collections import defaultdict
from CGAT import Experiment as E
from CGAT import GTF
from CGAT import IOTools
from CGAT.IndexedFasta import IndexedFasta
import textwrap
import pysam
import vcf
import re

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-b", "--bamfile", dest="bam", type="string",
                      help="BAM formated alignment file to test. Should have MD and NH tags set")
    parser.add_option("-t", "--quality-threshold", dest="threshold", type="int",
                       default=30,
                       help="minimum quality threshold for a mismatched base to count")
    parser.add_option("-f", "--fasta-path", dest="fastapath", type="string",
                       help="path to indexed fasta file for genome of choice")
    parser.add_option("-p", "--vcf-path", dest="vcfpath", type="string",
                       help="path to indexed vcf file for dataset  of choice")
    parser.add_option("-d", "--sample", dest="samppattern", type="string",
                       help="pattern to match and extract the donor name from the bam file, for use in parsing the vcf file")       
    #parser.add_option("-p", "--vcf-available", dest="vcfavail", type="int",
                       #help="whether a vcf file is available for this set of data")
    #parser.add_option("-a", "--vcf-avail", dest="vcfavail", type="int",
                       #help="integer indicating whether the VCF file provided was created by variant calling during the pipline (0) or VCF was already present (1)")
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    bamfile = pysam.AlignmentFile(options.bam)
    fastafile = IndexedFasta(options.fastapath)
    vcffile = vcf.Reader(open(options.vcfpath,"r"))
    options.stdout.write("\t".join(["gene_id",
                                    "mismatches",
                                    "bases",
                                    "low_qual",
                                    "a","t","c","g",
                                    "a_to_t","a_to_g","a_to_c",
                                    "t_to_a","t_to_g","t_to_c",
                                    "g_to_a","g_to_t","g_to_c",
                                    "c_to_a","c_to_t","c_to_g"]) + "\n")
    
    samplepattern = options.samppattern
    #samplepattern = samplepattern.replace('"','')
    #filed = open("test.txt","w")
    #filed.write("%s"%(samplepattern))
    #filed.close
    donorfrombam = re.search(r"%s"%(samplepattern),options.bam,flags=0).group(1)
#re.search("(\w+), (\w+)",str(deez[0]),flags=0).group(2)
    for gene in GTF.flat_gene_iterator(GTF.iterator(options.stdin)):
        
        exontuple=GTF.asRanges(gene, "exon")
        start = min(e.start for e in gene)
        end = max(e.end for e in gene)  

        seq = fastafile.getSequence(gene[0].contig, "+", start, end)
        thischr = gene[0].contig.replace("chr","")
        reads = bamfile.fetch(gene[0].contig, start, end)
        #filed = open("test.txt","w")
        #filed.write("%s\n"%(thischr)) 
        #filed.close      
        if all("chr" in c for c in vcffile.contigs.keys()) == False:            
            contig = (gene[0].contig).replace("chr","")
            if contig == "M":
                contig = contig + "T"
        else:
            contig = gene[0].contig
        #filed = open("test.txt","w")
        #filed.write("%s\n"%(contig))
        #filed.close       
        
        vcfregion = vcffile.fetch(contig,start,end)
        
        donors=[]
        regionchecker=[]
        for regionsnp in vcfregion:
            regionchecker.append(regionsnp)
        if regionchecker != []:
            samples = regionchecker[0].samples
            for dnr in samples:
                donors.append(dnr.sample)
                for samp in donors:
                    if (donorfrombam in samp) == True:
                        donorid = samp
    
                         
        gene_id = gene[0].gene_id
        mm_count = 0
        base_count = 0
        skipped = 0
        indelcount = 0
        matched_bases = defaultdict(int)
        transition = {"a_to_t":0,"a_to_g":0,"a_to_c":0,"t_to_a":0,"t_to_g":0,
        "t_to_c":0,"g_to_a":0,"g_to_t":0,"g_to_c":0,"c_to_a":0,"c_to_t":0,
        "c_to_g":0}
        
        snp_dict={}
        if regionchecker != []:
            for snp in vcfregion:
                if snp.genotype(donorid)["GT"] != "0/0":
                    snp_dict[snp.POS] = snp.ALT
        stepper=0    
        for read in reads:

            if read.is_unmapped:
                continue
            if read.is_duplicate:
                continue
            if read.mate_is_unmapped:
                continue
            if read.get_tag("NH") > 1:
                continue
            qualities = read.query_qualities

            alignmentcigar = read.cigarstring

            indelcount += (alignmentcigar.count("I") + alignmentcigar.count("D"))

            alignment = read.get_aligned_pairs(with_seq=True)

            testalignment = read.get_aligned_pairs(with_seq=True)
            
            def _is_exon_range(base):
                result_ranges=[]      
                for exonrange in exontuple:
                    result_ranges.append(exonrange[0] <= base[1] < exonrange[1]) 
                    if True in result_ranges:
                        return True
                    else:
                        return False
              

            alignment = [base for base in alignment 
                         if not base[0] is None and 
                         not base[1] is None
                         and _is_exon_range(base)]

           
            base_count += sum(1 for base in alignment
                          if start <= base[1] < end and
                          base[2].lower() != "n")
                          
	    #matched = [base for base in alignment
 		       #if not base[2].islower() and
                       #start <= base[1] < end]            
          
            #for base in matched:
                #if seq[(base[1])-start].lower() != base[2].lower():
                       #print read.query_alignment_sequence
                       #print seq[(alignment[0][1]-start):(alignment[-1][1]-start)]
                       #print seq[((base[1]-10)-start):((base[1]+10)-start)].lower()
                       #print read.tostring(bamfile)
                       #print start, end
                       #print seq[(base[1])-start]
                       #print base[2]
                       #print base[0]
                       #print base[1]-start                      
                       #print seq[(alignment[0][1]-start):(alignment[-1][1]-start)].upper()[base[0]]
                       #print ((alignment[0][1]-start) + base[0])
                       #print seq[(base[1])-start]
                       #print alignment[0][1]-start
                       #print read.get_aligned_pairs(with_seq=True)
                       #print textwrap.fill(seq,50)
                       #raise ValueError
                #else:
		       #matched_bases[base[2].lower()] += 1 
            
            total_alignment = [base for base in alignment
                               if start <= base[1] < end and
                               base[2].lower() != "n"]

            #for base in total_alignment:
                #if seq[(base[1])-start].lower() != base[2].lower():
                       #print read.query_alignment_sequence
                       #print seq[(alignment[0][1]-start):(alignment[-1][1]-start)]
                       #print seq[((base[1]-10)-start):((base[1]+10)-start)].lower()
                       #print read.tostring(bamfile)
                       #print start, end
                       #print seq[(base[1])-start]
                       #print base[2]
                       #print base[0]
                       #print base[1]-start                      
                       #print seq[(alignment[0][1]-start):(alignment[-1][1]-start)].upper()[base[0]]
                       #print ((alignment[0][1]-start) + base[0])
                       #print seq[(base[1])-start]
                       #print alignment[0][1]-start
                       #print read.get_aligned_pairs(with_seq=True)
                       #print textwrap.fill(seq,50)
                       #raise ValueError
                           
                #else:
                    #matched_bases[base[2].lower()] += 1

            for base in total_alignment:
                if seq[(base[1])-start].lower() != base[2].lower():
                    if (testalignment[0][1] is None) or (testalignment[-1][1] is None):
                        E.debug("first or last base of read is None")
                        raise ValueError
                    else:    
                        E.debug("identity of error causing base from read sequence: %s" %(read.query_alignment_sequence)[base[0]].lower())
                        E.debug("read sequence: %s" %(read.query_alignment_sequence))
                        E.debug("identity start and end of read as calculated from start and end as described in gtffile and extracted from fasta: %s" %(seq[(testalignment[0][1]-start):(testalignment[-1][1]-start)]))
                        E.debug("section of the read 10 bp downstream and upstream of the sequence containing the error extracted from the fasta: %s" %(seq[((base[1]-10)-start):((base[1]+10)-start)].lower()))
                        E.debug("filename?: %s" %(read.tostring(bamfile)))
                        E.debug("positions of start and end of the gene based on the gtf: %s,%s" %(start, end))
                        E.debug("identity of start of gene extratced from gtf: %s" %(seq[(base[1])-start]))
                        E.debug("identity of error causing base from reference genome: %s" %base[2])
                        E.debug("position of base in read: %s" %base[0])
                        E.debug("position of base in genome: %s" %base[1])
                        E.debug("position of base in read as calculated from position of base in genome and and start from gtf: %s" %(base[1]-start))                      
                        E.debug("identity of error causing base (reference), calculated from fasta and testalignment info: %s" %(seq[(testalignment[0][1]-start):(testalignment[-1][1]-start)].upper()[base[0]]))
                        E.debug("position of base in read from first alignment genome base minus start plus position of base in in read, should equal position of base in read: %s" %((testalignment[0][1]-start) + base[0]))
                        E.debug("identity of error causing base (reference), calculated from fasta and position of base in genome from aligned pairs: %s" %(seq[(base[1])-start]))
                        E.debug("position of start base in genome from the alignment minus position of start base in genome from the gtf, should be zero: %s" %(alignment[0][1]-start))
                        E.debug("complete aligned pairs, unfiltered: %s" %(testalignment))
                        E.debug("full fasta sequence of read: %s" %(textwrap.fill(seq,50)))
                        raise ValueError
                           
                else:
                    matched_bases[base[2].lower()] += 1
            if read.get_tag("NM") == 0:
                continue
         
            # mismatches
            

            readseq = read.query_alignment_sequence

	    def _is_snp(base):
		if len(snp_dict.keys()) > 0:
                    if snp_dict.has_key(base[1]):
                        if readseq[base[0]] != snp.dict[base[1]]:
                             return True
                        else:
                             return False
                    else:
                        return True
                else:
                    return True
            
            def _is_indel(base):
                if (len(readseq) >= (base[0] + 5)):
                    if (len(seq) < (((base[1])-start) + 5)):
                        upperrange = len(seq)-(base[1]-start)
                        lowerrange = 5 - upperrange                                            
                        readindelwindow=readseq[(base[0] - lowerrange):(base[0] + upperrange)]
                        seqindelwindow=seq[(((base[1])-start) - lowerrange):(((base[1])-start)+ upperrange)]                
                        matchwindows=[]
                        for i in range(len(readindelwindow)):
                            try:
                                matchwindows.append((readindelwindow[i].lower()==seqindelwindow[i].lower()))
                            except IndexError:
                                print i
                                print readindelwindow
                                print seqindelwindow
                                print start
                                print lowerrange
                                print upperrange
                                print base[0]
                                print (base[0] - lowerrange)
                                print (base[0] + upperrange)
                                print base[1]
                                print (base[1] - start)
                                print ((((base[1])-start) - lowerrange)-1)
                                print ((((base[1])-start) + upperrange)-1)
                                print readseq
                                print seq
                                print gene_id
                                print gene[0].contig                            
                                raise
                    elif (len(seq) >= (((base[1])-start) + 5)):
                        readindelwindow=readseq[base[0]:(base[0] + 5)]
                        seqindelwindow=seq[((base[1])-start):(((base[1])-start)+5)]
                        matchwindows=[]
                        for i in range(len(readindelwindow)):
                            try: 
                                matchwindows.append(readindelwindow[i].lower()==seqindelwindow[i].lower())
                            except IndexError:
                                print i
                                print readindelwindow
                                print seqindelwindow
                                print start
                                print base[0]
                                print base[1]
                                print (base[1] - start) - 1
                                print ((base[1] - start) + 5) - 1
                                print readseq
                                print seq
                                print gene_id
                                print gene[0].contig
                                raise
                    if matchwindows.count(False) >= 4:
                        return False
                    else:
                        return True
                elif (len(readseq) < (base[0] + 5)):
                    if len(seq) < (((base[1])-start) + 5):
                        readsequpperrange = len(readseq)-base[0]
                        readseqlowerrange = 5 - readsequpperrange
                        sequpperrange = len(seq) - (base[1] - start)
                        seqlowerrange = 5 - sequpperrange
                        if readsequpperrange < sequpperrange:
                            upperrange = readsequpperrange
                            lowerrange = readseqlowerrange
                        elif sequpperrange < readsequpperrange:
                            upperrange = sequpperrange
                            lowerrange = seqlowerrange
                        elif sequpperrange == readsequpperrange:
                            upperrange = sequpperrange
                            lowerrange = seqlowerrange                        
                    elif ((base[1] - start) - 4) < 0:
                        return True
                    else:
                        upperrange = len(readseq)-base[0]
                        lowerrange = 5 - upperrange
                    readindelwindow=readseq[(base[0] - lowerrange):(base[0] + upperrange)]
                    seqindelwindow=seq[(((base[1])-start) - lowerrange):(((base[1])-start)+ upperrange)]                
                    matchwindows=[]
                    for i in range(len(readindelwindow)):
                        try:
                            matchwindows.append((readindelwindow[i].lower()==seqindelwindow[i].lower()))
                        except IndexError:
                            print i
                            print readindelwindow
                            print seqindelwindow
                            print start
                            print lowerrange
                            print upperrange
                            print base[0]
                            print (base[0] - lowerrange)
                            print (base[0] + upperrange)
                            print base[1]
                            print (base[1] - start)
                            print ((((base[1])-start) - lowerrange))
                            print ((((base[1])-start) + upperrange))
                            print readseq
                            print seq
                            print gene_id
                            print gene[0].contig
                            
                            raise
                    if matchwindows.count(False) >= 4:
                        return False
                    else:
                        return True

            mismatches = [base for base in alignment
                          if base[2].islower() and
                          start <= base[1] < end and
		          qualities[base[0]] >= options.threshold and
                          base[2].lower() != "n" and 
                          _is_snp(base) and
                          readseq[base[0]].lower() != "n"
                          and _is_indel(base)]

            
            total_mm = sum(1 for base in alignment
                        if base[2].islower() and
                        start <= base[1] < end and
                        base[2].lower() != "n" and 
                        _is_snp(base) and
                        readseq[base[0]].lower() != "n")
            
            hq_mm = sum(1 for base in mismatches
                        if qualities[base[0]] >= options.threshold and
                        base[2].lower() != "n")


            for base in mismatches:
                genomebase = base[2].lower()
                readbase = readseq[base[0]].lower()
                try:
                    transition["%s_to_%s"%(genomebase, readbase)] += 1
                except KeyError:
                    print transition
                    print read.query_alignment_sequence.upper() 
                    print seq[(alignment[0][1]-start):(alignment[-1][1]-start)].upper()
	            print read.tostring(bamfile)
                    raise

	    
		    
            mm_count += hq_mm
            skipped += total_mm - hq_mm

        outline = "\t".join(map(str,[gene_id,
                                     mm_count,
                                     base_count,
                                     skipped,
                                     matched_bases['a'],
                                     matched_bases['t'],
                                     matched_bases['c'],
                                     matched_bases['g'],
                                     transition['a_to_t'],
                                     transition['a_to_g'],
                                     transition['a_to_c'],
                                     transition['t_to_a'],
                                     transition['t_to_g'],
                                     transition['t_to_c'],
                                     transition['g_to_a'],
                                     transition['g_to_t'],
                                     transition['g_to_c'],
                                     transition['c_to_a'],
                                     transition['c_to_t'],
                                     transition['c_to_g'],
                                     indelcount]))
        options.stdout.write(outline + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
