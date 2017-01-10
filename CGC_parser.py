#!/usr/bin/env python

################################################################
#
# CGC_parser.py
#
# Programmer:  Carol Zhou
#
# Description:  This code inputs the name of a gene caller plus
#    the gene-caller's output file, and outputs a properly formatted
#    file for input to CGC_main.py.  The purpose of this code is to
#    normalize all gene caller outputs to a standard format, which is
#    recognized by CGC_main.py (a code that compares gene calls
#    among gene callers.) Note that regardless of how the gene caller
#    designates the start/end of a gene call, this code sets the 
#    "leftEnd" to the smaller position number, and "rightEnd" to the
#    larger position number.  The strand ('+' or '-') determines 
#    whether these numbers represent start or end positions in the
#    script output file. 
#
# Input Files:  For GeneMarkS, use the XXX.fasta.lst file; for
#    Glimmer2, use the XXX.g2.coord file; and for Prodigal, use the
#    XXX.genes.sco file, for Glimmer3, use the run3.coords file. 
#
# Updates:
#    18 May 2016: begin
#    07 Jul 2016: Parses Glimmer3, Prodigal, GenemarkS, RAST
#    15 Aug 2016: upgraded to include PHATE parser
#
# Programmer's Notes:
#
################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.pdf FOR DETAILS.


import sys
import os
import re
import string
import copy
from subprocess import call

##### CONFIGURABLE

# GLIMMER3 bool controls which version of glimmer was run
# You get many more gene call matches when you adjust (by 3) the glimmer2 coordinates
# Glimmer3 changed this, so it should match better to other gene calls 
# You may adjust the Prodigal setting for the Prodigal out file you are using (sco vs. gff)
# RAST needs further testing; we are not running RAST for phage genomes

GLIMMER3 = True       # if False => glimmer2
PRODIGAL_sco = True   # using the XXX.genes.sco file
PRODIGAL_gff = False  # using the XXX.genes.gff file
RAST_GFF3 = True      # using RAST gff3 file; other RAST formats not yet supported

##### FILES

CODE_BASE = "./CGC_parser"
CODE      = CODE_BASE + ".py"
logfile   = CODE_BASE + ".log"
outfile   = CODE_BASE + ".out"
infile    = ""  # user provided
USER_OUT_PROVIDED = False 

LOGFILE = open(logfile,"w")

##### PATTERNS

p_comment  = re.compile('^#')
p_prodigal = re.compile('[Pp][Rr][Oo][Dd][Ii][Gg][Aa][Ll]')
p_genemark = re.compile('[Gg][Ee][Nn][Ee][Mm][Aa][Rr][Kk]')
p_glimmer  = re.compile('[Gg][Ll][Ii][Mm]+[Ee][Rr]')
p_rast     = re.compile('[Rr][Aa][Ss][Tt]')
p_phate    = re.compile('[Pp][Hh][Aa][Tt][Ee]')

##### IDIOMS

CHATTY = True
#CHATTY = False

##### CONSTANTS

HELP_STRING = "Script " + CODE + " inputs the name of a gene caller plus the output file arising \nfrom that gene caller. Then, the script converts the data to a format that is acceptable as input to \nscript CGC_main.py, which compares gene calls among a set of gene caller outputs.\nType: python" + CODE + " usage|input for more information\n"

USAGE_STRING = "Usage:  python " + CODE + " <geneCaller_name> <geneCall_filename> (optional)<output_filename>\n"

INPUT_STRING = "You may enter the name of a gene caller (e.g., Prodigal, GeneMark, Glimmer, RAST, PHATE), followed by the gene-call file that the program produced. For Prodigal, use the Name.genes.sco file. For GeneMarkS, use the Name.fasta.lst file. For Glimmer2, use the Name.g2.coord file, but for Glimmer3 use the run3.coords file. For RAST, use gff3 output. For PhATE... TBD.\n"

ACCEPTABLE_ARG_COUNT = (2,3,4)  # 2 if 'help'|'usage'|'input', or 3 if gene-caller and gene-caller.out, 4 if optional output file

##### GET INPUT PARAMETERS

geneCaller    = ""
geneCallerOut = ""
userOutfile   = ""

argCount = len(sys.argv)
if argCount in ACCEPTABLE_ARG_COUNT:
    match = re.search("help", sys.argv[1].lower())
    if match:
        print HELP_STRING
        LOGFILE.close(); exit(0)
    match = re.search("input", sys.argv[1].lower())
    if match:
        print INPUT_STRING
        LOGFILE.close(); exit(0)
    match = re.search("usage", sys.argv[1].lower())
    if match:
        print USAGE_STRING
        LOGFILE.close(); exit(0)

    # Capture name of gene caller and its output file
    if argCount == 3 or argCount == 4:
        geneCaller    = sys.argv[1].lower()  # case insensitive
        geneCallerOut = sys.argv[2]
    else:
        print USAGE_STRING
        LOGFILE.write("%s\n" % ("Incorrect number of command-line arguments provided"))
        LOGFILE.close(); exit(0)
    if argCount == 4:
        USER_OUT_PROVIDED = True
        userOutfile   = sys.argv[3]
else:
    print USAGE_STRING
    LOGFILE.write("%s\n" % ("Incorrect number of command-line arguments provided"))
    LOGFILE.close(); exit(0)

# Open file

fileError = False

try:
    INFILE = open(geneCallerOut,"r")
except IOError as e:
    fileError = True
    print e

try:
    OUTFILE = open(outfile,"w")
except IOError as e:
    fileError = True
    print e

if fileError:
    LOGFILE.write("%s%s%s\n" % ("ERROR: problem with input file:",geneCallerOut,e))
    LOGFILE.close(); exit(0)

if USER_OUT_PROVIDED:
    try:
        USER_OUT = open(userOutfile,"w")
    except IOError as e:
        fileError = True
        print e

##### FUNCTIONS

def ProcessGenemark(fLines,OUT):
    geneNo = 0; contig = ''; strand = ''; leftEnd = ''; rightEnd = ''; length = 0; count = 0; cclass = '' 
    p_dataLine   = re.compile('\s+(\d+)\s+([+-])\s+([\d\>\<]+)\s+(\d+)\s+(\d+)\s+\d+')
    p_contigLine = re.compile('FASTA\sdefinition\sline:\s(.*)\s+length=\d+')
    for line in fLines:
        match_comment    = re.search(p_comment,line)
        match_contigLine = re.search(p_contigLine,line)
        match_dataLine   = re.search(p_dataLine,line)
        if match_comment:
            continue 
        elif match_contigLine:
            contig = match_contigLine.group(1)
        elif match_dataLine:
            geneNo   = int(match_dataLine.group(1))
            strand   =     match_dataLine.group(2)     
            leftEnd  =     match_dataLine.group(3)   # Note: left/right end could have '<' or '>' symbol
            rightEnd =     match_dataLine.group(4)
            length   = int(match_dataLine.group(5))
            #cclass   =     match_dataLine.group(6)  # skip for now
            # Note: left/right end could have '>' or '<' symbol, if gene call spanned across contigs
            # I am removing the symbol, so if start=1 this may imply a gene that spans 
            # For a gene spanning 2 contigs or wrapping around, genemark estimates length based on
            # raw start/stop without symbol
            if re.search('^>',leftEnd) or re.search('^<',leftEnd):
                leftEnd = leftEnd[1:] 
            if re.search('^>',rightEnd) or re.search('^<',rightEnd):
                rightEnd = rightEnd[1:] 
            if strand != '+' and strand != '-':
                LOGFILE.write("%s%s\n" % ("ERROR: unknown strand designator, ",strand))
                OUT.write("%s\n" ("ERROR encountered: unknown strand designator\n"))
                if USER_OUT_PROVIDED:
                    USER_OUT.write("%s\n" ("ERROR encountered: unknown strand designator\n"))
                print "ERROR: unexpected strand designator,", strand
                return
            if contig == '':
                contig = 'unknown'  # Contig name may be absent in input file
            OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig))
            if USER_OUT_PROVIDED:
                USER_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig))
    return
            
def ProcessGlimmer(fLines,OUT):
    # NOTE:  Glimmer2 appears to truncate the initial met (sometimes?): test with other data sets (e.g., + strand calls)
    geneNo = 0; contig = ''; strand = ''; leftEnd = ''; rightEnd = ''; length = 0; count = 0; 
    p_contigLine = re.compile('^>(.*)\s+length=\d+\s+numreads=\d+') # Matches only for glimmer3

    if GLIMMER3:
        p_dataLine   = re.compile('orf(\d+)\s+(\d+)\s+(\d+)\s+([+-])\d\s+([\d\.]+)')
    else:
        p_dataLine = re.compile('\s+(\d+)\s+(\d+)\s+(\d+)\s+\[([+-])\d\sL=\s*(\d+)\sr=.*\]') 

    for line in fLines:
        match_comment    = re.search(p_comment,line)
        match_contigLine = re.search(p_contigLine,line)
        match_dataLine   = re.search(p_dataLine,line)
        if match_comment:
            continue 
        elif match_contigLine:
            contig = match_contigLine.group(1)
        elif match_dataLine:
            geneNo   = int(match_dataLine.group(1))
            left     = int(match_dataLine.group(2))
            right    = int(match_dataLine.group(3))
            strand   =     match_dataLine.group(4)
            if strand == '+':
                leftEnd  = left 
                if GLIMMER3:
                    rightEnd = right
                else:
                    rightEnd = right + 3
            elif strand == '-':  
                if GLIMMER3:
                    leftEnd = right
                else:
                    leftEnd  = right - 3    
                rightEnd = left  
            else:
                LOGFILE.write("%s%s\n" % ("ERROR: unknown strand designator, ",strand))
                OUT.write("%s\n" ("ERROR encountered: unknown strand designator\n"))
                if USER_OUT_PROVIDED:
                    USER_OUT.write("%s\n" ("ERROR encountered: unknown strand designator\n"))
                print "ERROR: unexpected strand designator,", strand
                return

            length = rightEnd - leftEnd + 1
            if contig == '':    # contig name may be left out of input file
                contig = 'unknown'
            OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,str(leftEnd),str(rightEnd),str(length),contig))

            if USER_OUT_PROVIDED:
                USER_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,str(leftEnd),str(rightEnd),str(length),contig))
    return

def ProcessRAST(fLines,OUT):
    geneNo = 0; contig = ""; strand = ''; leftEnd = ''; rightEnd = ''; length = 0; count = 0; 
    if RAST_GFF3:
        p_dataLine = re.compile('(.*)\t.*\tCDS\t(\d+)\t(\d+)\t\.\t([+-])\t\d\t.*')
        for line in fLines:
            match_comment = re.search(p_comment,line)
            match_dataLine = re.search(p_dataLine,line) 
            if match_comment:
                continue
            elif match_dataLine:
                count += 1
                geneNo = count 
                contig    =     match_dataLine.group(1)
                leftEnd   = int(match_dataLine.group(2))
                rightEnd  = int(match_dataLine.group(3))
                strand    =     match_dataLine.group(4) 
                length = rightEnd - leftEnd + 1
                OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,str(leftEnd),str(rightEnd),str(length),contig))
                if USER_OUT_PROVIDED:
                    USER_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,str(leftEnd),str(rightEnd),str(length),contig))
    else:
        pass # Not using other RAST format (for now)
    return
 
def ProcessProdigal(fLines,OUT):
    geneNo = 0; contig = "unknown"; strand = ''; leftEnd = ''; rightEnd = ''; length = 0; count = 0; 

    if PRODIGAL_sco:  # using XXX.genes.sco file
        p_dataLine = re.compile('^>(\d+)_(\d+)_(\d+)_([+-])')
        for line in fLines:
            match_comment = re.search(p_comment,line)
            match_dataLine = re.search(p_dataLine,line)
            if match_comment:
                continue 
            elif match_dataLine:
                geneNo   = match_dataLine.group(1)
                leftEnd  = int(match_dataLine.group(2))
                rightEnd = int(match_dataLine.group(3))
                strand   =     match_dataLine.group(4)
                length   = int(rightEnd) - int(leftEnd) + 1 
            OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig))
            if USER_OUT_PROVIDED:
                USER_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig))

    else: # using XXX.genes file
        p_dataLine = re.compile('(.*)\t.*\tCDS\t(\d+)\t(\d+)\t([\d\.]+)\t([+-])\t.\t(.*)') 
        for line in fLines:
            match_comment = re.search(p_comment,line)
            match_dataLine = re.search(p_dataLine,line)
            if match_comment:
                continue 
            elif match_dataLine:
                count += 1
                geneNo   = count 
                contig   =     match_dataLine.group(1)
                leftEnd  = int(match_dataLine.group(2))
                rightEnd = int(match_dataLine.group(3))
                strand   =     match_dataLine.group(5)
                length   = int(rightEnd) - int(leftEnd) + 1
            OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig))
            if USER_OUT_PROVIDED:
                USER_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig))
    return

def ProcessPhate(fLines,OUT):      # SDSU code
    p_dataLine = re.compile('^(\d+)\t(\d+)\t([+-])')
    geneNo = 0; contig = "unknown"; strand = ''; leftEnd = ''; rightEnd = ''; length = 0; count = 0 
    for line in fLines:
        match_comment  = re.search(p_comment,line)
        match_dataLine = re.search(p_dataLine,line)
        if match_comment:
            continue
        if match_dataLine:
            count += 1
            geneNo = count
            strand = match_dataLine.group(3)
            if strand == '+':
                leftEnd  = int(match_dataLine.group(1))
                rightEnd = int(match_dataLine.group(2))
            else:
                rightEnd = int(match_dataLine.group(1))
                leftEnd  = int(match_dataLine.group(2))
            length = abs(rightEnd - leftEnd) + 1
            OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig))
            if USER_OUT_PROVIDED:
                USER_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig))
    return

##### BEGIN MAIN 

# First, determine which gene caller was used

match_glimmer  = re.search(p_glimmer,geneCaller)
match_genemark = re.search(p_genemark,geneCaller)
match_prodigal = re.search(p_prodigal,geneCaller)
match_rast     = re.search(p_rast,geneCaller)
match_phate    = re.search(p_phate,geneCaller)

OUTFILE.write("%s%s%s%s%s\n" % ('# ',geneCaller, " gene calls",", taken from file ",geneCallerOut))
OUTFILE.write("%s\n" % ("Gene No.\tStrand\tLeftEnd\tRightEnd\tLength\tContig"))

if USER_OUT_PROVIDED:
    USER_OUT.write("%s%s%s%s%s\n" % ('# ',geneCaller, " gene calls",", taken from file ",geneCallerOut))
    USER_OUT.write("%s\n" % ("Gene No.\tStrand\tLeftEnd\tRightEnd\tLength\tContig"))

fileLines = INFILE.read().splitlines()

if match_genemark:
    ProcessGenemark(fileLines,OUTFILE)
elif match_glimmer:
    ProcessGlimmer(fileLines,OUTFILE)
elif match_prodigal:
    ProcessProdigal(fileLines,OUTFILE)
elif match_rast:
    ProcessRAST(fileLines,OUTFILE)
elif match_phate:
    ProcessPhate(fileLines,OUTFILE)
else:
    LOGFILE.write("%s%s\n" % ("ERROR: Cannot process unknown gene caller output file:",geneCaller))

if USER_OUT_PROVIDED:
    USER_OUT.write("%s\n" % ("# END"))
OUTFILE.write("%s\n" % ("# END"))

##### CLEAN UP

INFILE.close()
OUTFILE.close()
if USER_OUT_PROVIDED:
    USER_OUT.close()
LOGFILE.write("%s\n" % ("Processing complete"))
LOGFILE.close()
