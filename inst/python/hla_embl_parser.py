# -*- coding: cp1252 -*-
'''
Created on 11.09.2015

parser functions for repeat use

@author: schoene
'''
#===========================================================
#import modules:

##from sys import argv
from pickle import dump
import RPython

#===========================================================
#global definitions
#global is_integer
#global read_dat_file

#===========================================================
#classes:

class Allele:
    def __init__(self, ID, locus, name, seq, length, UTR5, UTR3, exon_dic, intron_dic, exonpos_dic, intronpos_dic, utrpos_dic, target):
        self.ID = ID # Accession number or name
        self.locus = locus # gene
        self.name = name # allele name
        self.seq = seq.upper() # full sequence in upper case
        self.length = len(seq) # length of sequence
        self.UTR3 = UTR3 # sequence of UTR3
        self.UTR5 = UTR5 # sequence of UTR5
        self.exon_dic = exon_dic # dict of format {1:'EXON1_SEQ', 2:'EXON2_SEQ'...}
        self.intron_dic = intron_dic # dict of format {1:'INTRON1_SEQ', 2:'INTRON2_SEQ'...}
        self.exonpos_dic = exonpos_dic
        self.intronpos_dic = intronpos_dic
        self.utrpos_dic = utrpos_dic 
        self.is_ref = False
        self.full_seq = False # True if at least some introns are known, otherwise false
        if len(intron_dic.keys()) > 0:
            self.full_seq = True
        
        # if name.startswith("HLA"):
        #     name = name.split("-")[1]
        self.fasta_header = ">%s" % name
        #self.fasta_header = ">%s %s" % (name, self.length)
        # calculate CDS:
        exons = exon_dic.keys()
        exons.sort()
        self.CDS = ""
        for exon in exons:
            self.CDS += exon_dic[exon]
        
    def __repr__(self):
        return self.name

 
#===========================================================
# reading functions:

def read_dat_file(dat_file, target, isENA = False, verbose = False):
    """reads content of a .dat file (EMBL format),
    returns list of allele objects.
    The parameter 'target' expects one of the following: "HLA", "Blutgruppen","CCR5", "KIR".
    """
    alleles = []
    if verbose:
        print "Lese %s ein..." % dat_file
    with open(dat_file, "r") as f:
        data = f.readlines()
        for i in range(len(data)):
            line = data[i]
            if line.startswith("ID"):
                s = line.split()
                allele_ID = s[1].replace(";","")
                length = s[-2] ##was 7 in orignial code but it seems to be constitently the second last (stefan: 2016_10_05)
                seq = ""
                UTR3 = ""
                UTR5 = ""
                UTR5_start = False
                UTR5_end = False
                UTR3_start = False
                UTR3_end = False
                exon_dic = {}
                exonpos_dic = {}
                intron_dic = {}
                intronpos_dic = {}
                utrpos_dic = {} 
                
                if target in ["Blutgruppen","CCR5"]:
                    allele = allele_ID
                    if allele.find("*")>0:
                        locus = allele.split("*")[0]
                    elif allele.find("_")>0:
                        locus = allele.split("_")[0]
                    else:
                        print "!!!Cannot see Locus of Allele %s! Please adjust Input file!" % allele
                        print line
                        sys.exit()
            
            elif line.startswith("DE"):
                s = line.split()
                if target in ["HLA", "HLA_23_with_introns", "Phasing_HLA_23", "KIR"]:
                    allele = s[1][:-1]
                    locus = allele.split("*")[0]

            elif line.startswith("FT"):
                s = line.split()
                if s[1] == "UTR":
                    start = int(s[-1].split(".")[0]) - 1
                    end = int(s[-1].split(".")[-1])

                    if start == 0:
                        UTR5_start = start
                        UTR5_end = end
                        utrpos_dic["utr5"] = (start, end)
                    else:
                        UTR3_start = start
                        UTR3_end = end
                        utrpos_dic["utr3"] = (start, end)
                elif s[1] == "exon":
                    start = int(s[-1].split(".")[0]) - 1
                    end = int(s[-1].split(".")[-1])
                    next_line = data[i+1]
                    assert next_line.find("number") > 0, "Cannot find exon number in %s:\n '%s'\n '%s'" % (allele, line, next_line)
                    exon_num = int(next_line.split('"')[-2])
                    exonpos_dic[exon_num] = (start, end)

                elif s[1] == "intron":
                    start = int(s[-1].split(".")[0]) - 1
                    end = int(s[-1].split(".")[-1])
                    next_line = data[i+1]
                    assert next_line.find("number") > 0, "Cannot find intron number in %s:\n '%s'\n '%s'" % (allele, line, next_line)
                    intron_num = int(next_line.split('"')[-2])
                    intronpos_dic[intron_num] = (start, end)

            elif line.startswith("SQ"):
                read_on = True
                j = 0
                while read_on:
                    j += 1
                    s = data[i+j]
                    if s.startswith("//"):
                        read_on = False
                    else:
                        if is_integer1(s.split()[-1]):
                            myseq = "".join(s.split()[:-1]).upper()
                        else:
                            myseq = "".join(s.split()).upper()
                        seq += myseq

            elif line.startswith("//"):
                for exon in exonpos_dic:
                    (start,end) = exonpos_dic[exon]
                    exon_seq = seq[start:end]
                    exon_dic[exon] = exon_seq

                for intron in intronpos_dic:
                    (start,end) = intronpos_dic[intron]
                    intron_seq = seq[start:end]
                    intron_dic[intron] = intron_seq

                if UTR5_end:
                    UTR5 = seq[UTR5_start:UTR5_end].upper()
                    #print(UTR5)
                if UTR3_end:
                    UTR3 = seq[UTR3_start:UTR3_end].upper()
                    #print(UTR3)

                myAllele = Allele(allele_ID, locus, allele, seq, length, UTR5, UTR3, exon_dic, intron_dic, exonpos_dic, intronpos_dic, utrpos_dic, target)
                if target == "HLA": # HLA.dat contains other loci, too - MIC, TAP...
                    usable = False
                    usable_loci = ["HLA-A*", "HLA-B*", "HLA-C*", "HLA-DPB1*", "HLA-DQB1*", "HLA-DRB"]
                    for loc in usable_loci:
                        if allele.startswith(loc): # HLA.dat contains other loci, too - MIC, TAP...
                            usable = True
                    if usable:
                        alleles.append(myAllele)
                else:
                    alleles.append(myAllele)
    if verbose:
        print "\t%s Allele von %s erfolgreich eingelesen!" % (len(alleles), target)

    alleleHash = {}
    for allele in alleles:
	if (not isENA and (allele.name.find("DQB1") != -1)) and \
	   ((allele.name != "HLA-DQB1*05:03:01:01") or (allele.name != "HLA-DQB1*06:01:01")): continue
	alleleHash[allele.name] = allele
    

    return alleleHash


def read_dat_file_simple(dat_file, target, isENA = False, verbose = False):
    returnVec = []
    returnName= []
    myHash = read_dat_file(dat_file, target, isENA, verbose)
    for k in myHash.keys():
        h = myHash[k]
        returnVec.append(getSimple(h))
        returnName.append(h.locus)

    returnVec.append(returnName)
    return returnVec


def read_dat_file_simple_locus(dat_file, target, gene, isENA = False, verbose = False):
    returnVec = []
    returnName= []
    myHash = read_dat_file(dat_file, target, isENA, verbose)
    for k in myHash.keys():
        if (myHash[k].locus==gene):
            h = myHash[k]
            returnVec.append(getSimple(h))
            returnName.append(h.locus)

    returnVec.append(returnName)
    return returnVec




def is_integer1(x):
    try:
        int(x)
        return True
    except ValueError:
        return False


def getSimple(myAllele):
    ##def getSimple(myAllele,myReturn):
    myReturn = [ ]
    ##myReturn.append("NEW")
    myReturn.append(myAllele.ID)
    myReturn.append(myAllele.locus)
    myReturn.append(myAllele.name)
    myReturn.append(myAllele.seq)
    myReturn.append(myAllele.length)
    myReturn.append(myAllele.is_ref)
    myReturn.append(myAllele.full_seq)
    for utr in myAllele.utrpos_dic.keys():
        (start,end) = myAllele.utrpos_dic[utr]
        myReturn.append(utr+" "+str(start+1)+" "+str(end))
        ##note for R use indexes starting with 1 (1-based), and the end is inclusive
    for exon in myAllele.exonpos_dic.keys():
        (start,end) = myAllele.exonpos_dic[exon]
        myReturn.append("Exon "+str(start+1)+" "+str(end)+" "+str(exon))
        ##note for R use indexes starting with 1 (1-based), and the end is inclusive
    for intron in myAllele.intronpos_dic.keys():
        (start,end) = myAllele.intronpos_dic[intron]
        myReturn.append("Intron "+str(start+1)+" "+str(end) + " "+str(intron))
        ##note for R use indexes starting with 1 (1-based), and the end is inclusive    
    
    return RPython.vectorR(myReturn, "character")
    ##return(myReturn)

#===========================================================
# writing functions:

def write_fasta(alleles, output_fasta, verbose = False):
    """takes a list of allele objects,
    writes a fasta-file containing their full sequences
    """
    if verbose:
        print "Schreibe %s Allele nach %s..." % (len(alleles), output_fasta)
    with open(output_fasta, "w") as g:
        for k in alleles.keys():
            g.write("%s\n%s\n" % (alleles[k].fasta_header, alleles[k].seq))
    if verbose:
        print "\tFertig!"

if __name__ == '__main__':

    imgt_download_file, target, outfile = argv[1], argv[2], argv[3]
    alleles = read_dat_file(imgt_download_file, target)
   
    alleleNames = alleles.keys()
    print "Total Alleles: ", len(alleleNames)
   
    counts = {"utr3":0,"utr5":0,"both":0}
    for alleleName in alleleNames:
	alleleObj = alleles[alleleName]
	if (len(alleleObj.UTR3) & len(alleleObj.UTR5)): counts["both"] += 1
	elif (len(alleleObj.UTR3)): counts["utr3"] += 1
	elif (len(alleleObj.UTR5)): counts["utr5"] += 1
	else: pass
	
    print "full length : ", counts["both"]
    print "utr3 only : ", counts["utr3"]
    print "utr5 only : ", counts["utr5"]

    write_fasta(alleles, outfile)

   #alleleDumpFile = open(argv[3], "w")
   #dump(alleles, alleleDumpFile)
   #alleleDumpFile.close()  
