#!/usr/bin/env python

################################################################
#
# CGC_main.py  # Compare Gene Calls Main
#
# Programmer: Carol Zhou
#
# Description:  Accepts a list of input files, each comprising a set of
#    gene calls from a given gene caller program (e.g., Prodigal, Glimmer,
#    GeneMark, PhATE).  Outputs comparisons accross the gene calls. 
#    Note:  The input files are re-formatted using CGC_parser.py, so that
#    they have a common format and identify the gene caller in the comments
#    at the top of the file.
#
# Updates:
#    17 May 2016: begin
#    3 June 2016: adding CGC_geneCall.py
#    21 June 2016: ready for code release
#    6 July 2016: added RAST; added Prodigal gff; capturing contig when given
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
import CGC_geneCall
import CGC_compare

##### FILES

CODE_BASE = "./CGC_main"
CODE_FILE = CODE_BASE + ".py"
LOG_FILE  = CODE_BASE + ".log"
OUT_FILE  = CODE_BASE + ".out"

infile = ""
OUT = open(OUT_FILE,"w") 
LOG = open(LOG_FILE,"w")

##### PATTERNS

p_comment  = re.compile('^#')
p_order    = re.compile('Order')

##### PRINT CONTROL 

CHATTY = True  # This will print status as the code executes
#CHATTY = False

#DEBUG = True    # Print even more!
DEBUG = False

##### CONSTANTS

HELP_STRING = "This code inputs a list of at least 2 files comprising gene calls (generated by a gene caller program) and outputs the genes that are in common and unique with respect to each caller.  Type: python " + CODE_FILE + " usage|input|detail for more information\n"

USAGE_STRING = "Usage:  python " + CODE_FILE + " <infile>n\n"

INPUT_STRING = "Input for " + CODE_FILE + " comprises a list of path/filenames comprising outputs generated by gene caller programs, with each separated by a single space. The output files are to have been prepared using code \"CGC_parser.py\" to assure that they have a common pre-defined format and indicate the name of the gene caller in the comments section.\nExample:  python " + CODE_FILE + " genemark.calls prodigal.calls\n"

INFO_STRING = "This code currently supports the following gene callers:  GeneMark, Glimmer, Prodigal, RAST, and PhATE. For more information regarding input to " + CODE_FILE + ", type:  " + CODE_FILE + " input"

##### GET INPUT PARAMETERS

fileSet = []
argCount = len(sys.argv)
if argCount > 1:
    match = re.search("help", sys.argv[1].lower())
    if match:
        print HELP_STRING
        LOG.close(); exit(0)
    match = re.search("input", sys.argv[1].lower())
    if match:
        print INPUT_STRING
        LOG.close(); exit(0)
    match = re.search("usage", sys.argv[1].lower())
    if match:
        print USAGE_STRING
        LOG.close(); exit(0)
    match = re.search("detail", sys.argv[1].lower())
    if match:
        print INFO_STRING
        LOG.close(); exit(0)
    match = re.search("info", sys.argv[1].lower())
    if match:
        print INFO_STRING
        LOG.close(); exit(0)
    else:
        fileSet = sys.argv[1:]  # skip 0th element = name of code
else:
    LOG.write("%s\n" % ("Incorrect number of command-line arguments provided"))
    print USAGE_STRING
    LOG.close(); exit(0)

##### BEGIN MAIN 

count = 0
callerList = []
callSet_obj = CGC_geneCall.GeneCallSet()

# For each user-provided gene call file, create a call set and add to list of call sets

if CHATTY:
    print "Main: Iterating through fileSet..."

for geneFile in fileSet:
    callSet = copy.deepcopy(callSet_obj)
    geneFile_handle = open(geneFile,"r")
    if CHATTY:
        print "Adding Calls from file", geneFile
    callSet.AddGeneCalls(geneFile_handle)
    geneFile_handle.close()
    callerList.append(callSet)

if CHATTY:
    print "Main: callerList is", 
    for caller in callerList:
        print caller.geneCaller, ', ',
    print 

# Check
if DEBUG:
    print "\n******************Original Lists:"
    for caller in callerList:
        caller.PrintAll()
        print

# Sort calls in each list

if CHATTY:
    print "Main: Sorting gene calls for each caller..."

for caller in callerList:
    caller.SortGeneCalls()

# Check
if DEBUG:
    print "\n******************Sorted Lists:"
    for caller in callerList:
        caller.PrintAll()
        print

# Compare across the call sets

if CHATTY:
    print "Main: Comparing accross the call sets..."

compareGCs = CGC_compare.Comparison()
for caller in callerList:
    compareGCs.Merge(caller.geneCallList)
compareGCs.Compare()
compareGCs.IdentifyCommonCore()
compareGCs.PrintReport()

# Check
if DEBUG:
    compareGCs.PrintAll()

##### CLEAN UP

OUT.close()
LOG.close()
