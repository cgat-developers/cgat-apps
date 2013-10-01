################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: PipelineGO.py 2877 2010-03-27 17:42:26Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""
======================================================
PipelineMappingQC.py - common tasks for QC'ing mapping
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------


Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----


"""

import CGAT.Experiment as E
import logging as L
import CGAT.Database as Database
import CGAT.CSV as CSV

import sys
import os
import re
import shutil
import itertools
import math
import glob
import time
import gzip
import collections
import random

import numpy
import sqlite3
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
import rpy2.rinterface as ri

import CGAT.Pipeline as P
import pysam

try:
    PARAMS = P.getParameters()
except IOError:
    pass

def getPicardOptions():
    return "-pe dedicated 3 -R y -l mem_free=1.4G -l picard=1"

def getNumReadsFromReadsFile( infile ):
    '''get number of reads from a .nreads file.'''
    with IOTools.openFile( infile ) as inf:
        line = inf.readline()
        if not line.startswith( "nreads" ):
            raise ValueError( "parsing error in file '%s': expected first line to start with 'nreads'")
        nreads = int( line[:-1].split("\t")[1] )
    return nreads

def getNumReadsFromBAMFile( infile ):
    '''count number of reads in bam file.'''
    # by-passes a problem with pysam, which was reading in stdout as the first elements in list data
    tmpf = P.getTempFile( ".")
    tmpfile_name = tmpf.name
    statement = '''samtools idxstats %(infile)s > %(tmpfile_name)s'''
    
    P.run()

    read_info = IOTools.openFile( tmpfile_name ).readlines()
    os.unlink( tmpfile_name )

    try:
        data = sum( map(int, [ x.split("\t")[2] for x in read_info if not x.startswith("#")]  ) )

    except IndexError, msg:
        raise IndexError( "can't get number of reads from bamfile, msg=%s, data=%s" % (msg, read_info))
    return data


def buildPicardInsertSizeStats( infile, outfile, genome_file ):
    '''gather BAM file insert size statistics using Picard '''

    to_cluster = True
    job_options = getPicardOptions()

    if getNumReadsFromBAMFile(infile) == 0:
        E.warn( "no reads in %s - no metrics" % infile )
        P.touch( outfile )
        return

    statement = '''CollectInsertSizeMetrics
                                       INPUT=%(infile)s 
                                       REFERENCE_SEQUENCE=%(genome_file)s
                                       ASSUME_SORTED=true 
                                       OUTPUT=%(outfile)s 
                                       VALIDATION_STRINGENCY=SILENT 
                   >& %(outfile)s'''

    P.run()

def buildPicardAlignmentStats( infile, outfile, genome_file ):
    '''gather BAM file alignment statistics using Picard '''

    to_cluster = True
    job_options = getPicardOptions()

    if getNumReadsFromBAMFile(infile) == 0:
        E.warn( "no reads in %s - no metrics" % infile )
        P.touch( outfile )
        return

    # Picard seems to have problem if quality information is missing
    # or there is no sequence/quality information within the bam file.
    # Thus, add it explicitely.
    statement = '''cat %(infile)s 
                       | python %(scriptsdir)s/bam2bam.py -v 0 --set-sequence --sam 
                       | CollectMultipleMetrics 
                                       INPUT=/dev/stdin 
                                       REFERENCE_SEQUENCE=%(genome_file)s
                                       ASSUME_SORTED=true 
                                       OUTPUT=%(outfile)s 
                                       VALIDATION_STRINGENCY=SILENT 
                       >& %(outfile)s'''

    P.run()

def buildPicardDuplicationStats( infile, outfile ):
    '''Record duplicate metrics using Picard, the marked records are discarded'''

    to_cluster = True
    job_options = getPicardOptions()

    if getNumReadsFromBAMFile(infile) == 0:
        E.warn( "no reads in %s - no metrics" % infile )
        P.touch( outfile )
        return

    statement = '''MarkDuplicates
                                   INPUT=%(infile)s 
                                   ASSUME_SORTED=true 
                                   METRICS_FILE=%(outfile)s
                                   OUTPUT=/dev/null 
                                   VALIDATION_STRINGENCY=SILENT 
                 '''
    P.run()

def buildPicardDuplicateStats( infile, outfile ):
    '''Record duplicate metrics using Picard and keep the dedupped .bam file'''

    to_cluster = True
    job_options = getPicardOptions()

    if getNumReadsFromBAMFile(infile) == 0:
        E.warn( "no reads in %s - no metrics" % infile )
        P.touch( outfile )
        return

    statement = '''MarkDuplicates
                                   INPUT=%(infile)s 
                                   ASSUME_SORTED=true 
                                   METRICS_FILE=%(outfile)s.duplicate_metrics
                                   OUTPUT=%(outfile)s 
                                   VALIDATION_STRINGENCY=SILENT 
                 '''
    P.run()

def buildPicardGCStats( infile, outfile, genome_file ):
    '''Gather BAM file GC bias stats using Picard '''
    to_cluster = True

    job_options = getPicardOptions()

    if getNumReadsFromBAMFile(infile) == 0:
        E.warn( "no reads in %s - no metrics" % infile )
        P.touch( outfile )
        return

    statement = '''CollectGcBiasMetrics
                                       INPUT=%(infile)s 
                                       REFERENCE_SEQUENCE=%(genome_file)s
                                       OUTPUT=%(outfile)s 
                                       VALIDATION_STRINGENCY=SILENT 
                                       CHART_OUTPUT=%(outfile)s.pdf 
                                       SUMMARY_OUTPUT=%(outfile)s.summary
                   >& %(outfile)s'''

    P.run()

def loadPicardMetrics( infiles, outfile, suffix, pipeline_suffix = ".picard_stats" ):
    '''load picard metrics.'''

    tablename = P.toTable( outfile )
    tname = "%s_%s" % (tablename, suffix)

    outf = P.getTempFile( "." )

    filenames = [ "%s.%s" % (x, suffix) for x in infiles ]

    first = True


    for filename in filenames:
        track = P.snip( os.path.basename(filename), "%s.%s" % (pipeline_suffix, suffix ) )

        if not os.path.exists( filename ): 
            E.warn( "File %s missing" % filename )
            continue

        lines = IOTools.openFile( filename, "r").readlines()
        
        # extract metrics part
        rx_start = re.compile("## METRICS CLASS")
        for n, line in enumerate(lines):
            if rx_start.search(line ):
                lines = lines[n+1:]
                break

        for n, line in enumerate(lines):
            if not line.strip(): 
                lines = lines[:n]
                break
            
        if len(lines) == 0:
            E.warn("no lines in %s: %s" % (track,f))
            continue
        if first: 
            outf.write( "%s\t%s" % ("track", lines[0] ) )
            fields = lines[0][:-1].split("\t")
        else:
            f = lines[0][:-1].split("\t")
            if f != fields:
                raise ValueError("file %s has different fields: expected %s, got %s" % \
                                     (filename, fields, f ))

        first = False
        for i in range(1, len(lines)):
            outf.write( "%s\t%s" % (track,lines[i] ))
            
    outf.close()

    tmpfilename = outf.name

    to_cluster = False
    statement = '''cat %(tmpfilename)s
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tname)s 
                      --allow-empty
                > %(outfile)s
               '''
    P.run()

    os.unlink( tmpfilename )

def loadPicardHistogram( infiles, outfile, suffix, column, pipeline_suffix = ".picard_stats" ):
    '''extract a histogram from a picard output file and load it into database.'''

    tablename = P.toTable( outfile )
    tname = "%s_%s" % (tablename, suffix)
    
    tname = P.snip( tname, "_metrics") + "_histogram"

    # some files might be missing
    xfiles = [ x for x in infiles if os.path.exists( "%s.%s" % (x, suffix) ) ]

    if len(xfiles) == 0: 
        E.warn ( "no files for %s" % tname )
        return
    
    header = ",".join( [P.snip( os.path.basename(x), pipeline_suffix) for x in xfiles ] )        
    filenames = " ".join( [ "%s.%s" % (x, suffix) for x in xfiles ] )

    # there might be a variable number of columns in the tables
    # only take the first ignoring the rest
    statement = """python %(scriptsdir)s/combine_tables.py
                      --regex-start="## HISTOGRAM"
                      --missing=0
                      --take=2
                   %(filenames)s
                | python %(scriptsdir)s/csv2db.py
                      --header=%(column)s,%(header)s
                      --replace-header
                      --index=track
                      --table=%(tname)s 
                >> %(outfile)s
                """
    
    P.run()

def loadPicardAlignmentStats( infiles, outfile ):
    '''load all output from Picard's CollectMultipleMetrics
    into sql database.'''

    loadPicardMetrics( infiles, outfile, "alignment_summary_metrics" )

    # insert size metrics only available for paired-ended data
    loadPicardMetrics( infiles, outfile, "insert_size_metrics" )

    histograms = ( ("quality_by_cycle_metrics", "cycle"),
                   ("quality_distribution_metrics", "quality"),
                   ("insert_size_metrics", "insert_size" ) )

    for suffix, column in histograms:
        loadPicardHistogram( infiles, outfile, suffix, column )


def loadPicardDuplicationStats( infiles, outfile ):
    '''load picard duplicate filtering stats.'''
    # SNS: added to enable naming consistency

    suffix = "duplication_metrics"

    # the loading functions expect "infile_name.pipeline_suffix" as the infile names.
    infiles = [ x[ 0 : -len( "." + suffix ) ] for x in infiles]
    
    loadPicardMetrics( infiles, outfile, suffix )
    loadPicardHistogram( infiles, outfile, suffix, "coverage_multiple" )

def loadPicardDuplicateStats( infiles, outfile ):
    '''load picard duplicate filtering stats.'''

    loadPicardMetrics( infiles, outfile, "duplicate_metrics", pipeline_suffix = ".bam" )
    loadPicardHistogram( infiles, outfile, "duplicate_metrics", "duplicates", pipeline_suffix = ".bam" )

def buildBAMStats( infile, outfile ):
    '''Count number of reads mapped, duplicates, etc. '''
    to_cluster = True

    statement = '''python %(scriptsdir)s/bam2stats.py 
                          --force 
                          --output-filename-pattern=%(outfile)s.%%s 
                          < %(infile)s 
                          > %(outfile)s'''
    P.run()


###########################################################
###########################################################
###########################################################
def loadBAMStats( infiles, outfile ):
    '''load bam2stats.py output into sqlite database.'''

    # scriptsdir = PARAMS["general_scriptsdir"]
    header = ",".join( [P.snip( os.path.basename(x), ".readstats") for x in infiles] )
    filenames = " ".join( [ "<( cut -f 1,2 < %s)" % x for x in infiles ] )
    tablename = P.toTable( outfile )
    E.info( "loading bam stats - summary" )
    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/track/"
                | perl -p -e "s/unique/unique_alignments/"
                | python %(scriptsdir)s/table2table.py --transpose
                | python %(scriptsdir)s/csv2db.py
                      --allow-empty
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s"""
    P.run()

    for suffix in ("nm", "nh"):
        E.info( "loading bam stats - %s" % suffix )
        filenames = " ".join( [ "%s.%s" % (x, suffix) for x in infiles ] )
        tname = "%s_%s" % (tablename, suffix)
        
        statement = """python %(scriptsdir)s/combine_tables.py
                      --header=%(header)s
                      --skip-titles
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/%(suffix)s/"
                | python %(scriptsdir)s/csv2db.py
                      --table=%(tname)s 
                      --allow-empty
                >> %(outfile)s """
        P.run()

    # load mapping qualities, there are two columns per row
    # 'all_reads' and 'filtered_reads'
    # Here, only filtered_reads are used (--take=3)
    for suffix in ("mapq",):
        E.info( "loading bam stats - %s" % suffix )
        filenames = " ".join( [ "%s.%s" % (x, suffix) for x in infiles ] )
        tname = "%s_%s" % (tablename, suffix)
        
        statement = """python %(scriptsdir)s/combine_tables.py
                      --header=%(header)s
                      --skip-titles
                      --missing=0
                      --ignore-empty
                      --take=3
                   %(filenames)s
                | perl -p -e "s/bin/%(suffix)s/"
                | python %(scriptsdir)s/csv2db.py
                      --table=%(tname)s 
                      --allow-empty
                >> %(outfile)s """
        P.run()

