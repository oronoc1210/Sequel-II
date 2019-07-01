#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
---------------------------------------------------------------------------------#
author: Jie Guo
creation date: 20180809
This is the script of iso-seq3 Pipeline

Version 2: 20190429
add log function and change run_command to run_cmd
'''

import os
import sys
import re
import platform
import time
import subprocess
import argparse
import pandas as pd
sys.path.insert(0, "/global/dna/projectdirs/PI/rqc/prod/jgi-rqc-pipeline/lib/")
from common import checkpoint_step, run_cmd, get_logger

'''
Function to check inputs and setup output dir
'''
def setup():
    
    global log
    log_file = 'iso-seq3.log'
    print_log = args.print_log
    log_level = "INFO"
    log = get_logger("peakcalling", log_file, log_level, print_log)

    log.info('setup')    
    #set global variables
    global program_execution
    global subreads
    global primers
    global n_threads
    global prefix
    global n_split
    global ccs_reads
    global reference
    global gtf
    global my_dir
    global samtools
    
    samtools = 'shifter --image=docker:rmonti/samtools samtools'
    ccs_reads = None
    reference = None
    gtf = None
    
    n_split = 36
    program_execution = ' '.join(sys.argv)    #cmd line
    subreads = os.path.abspath(args.subreads)    #read fasta
    primers = os.path.realpath(args.primers)    #reference/db fasta
    n_threads = args.num_threads
    prefix = args.prefix
    machine_info = ', '.join(platform.uname()).rstrip(', ') # machine platform info
    #n_cpu = multiprocessing.cpu_count() #machine cpu info
    my_dir = os.path.dirname(__file__)
    
    if args.ccs:
        ccs_reads = os.path.abspath(args.ccs)
    if args.reference:
        reference = os.path.realpath(args.reference)
    if args.gtf:
        if not args.gtf.endswith('gtf'):
            log.error('annotation file should be gtf format')
            sys.exit(1)
        else:
            gtf = os.path.realpath(args.gtf)
        
    #make output directory if not created
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)
        
    if not os.path.exists("%s.pbi" % (args.subreads)):
        #cmd = "shifter -i bryce911/smrtlink pbindex %s" % args.subreads
        cmd="shifter --image=registry.services.nersc.gov/jgi/smrtlink:5.0.1.9585 /smrtlink/smrtcmds/bin/pbindex %s" % args.subreads
        log.info('no pbi index, start to generate one')
        std_out, std_err, exit_code = run_cmd(cmd, log)
         #subprocess.call(cmd, shell=True)
    #chdir to output dir

    os.chdir(args.out_dir)
    
    #print slash info
    log.info('setup: %s: %s' % (my_name,version))
    log.info('setup: %s: %s' % ('prog_exec',program_execution))
    log.info('setup: %s: %s' % ('output_dir',args.out_dir))
    log.info('setup: %s: %s' % ("subreads", subreads))
    log.info('setup: %s: %s' % ("primers", primers))
    log.info('setup: %s: %s' % ("run_machine", machine_info))
    #log('setup',"%s: %s" % ("avail_cpu", n_cpu))


def run_CCS():
    log.info('start to run ccs on subreads')
    cmd = "cp -f %s/resolved-tool-contract.json ." % my_dir
    std_out, std_err, exit_code = run_cmd(cmd, log)
    ccs_reads = "%s.ccs.bam" % prefix
    script = "ccs"
    option = "--minReadScore=0.65 --minZScore=-3.5 --numThreads=32 --minSnr=4 --maxLength=21000 --minLength=10 --minPasses=3 --minPredictedAccuracy=0.98 --maxDropFraction=0.34" % (n_threads)
    cmd = "%s %s %s %s " % (script, option, subreads, ccs_reads)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    ccs_reads = os.path.realpath(ccs_reads)
    return ccs_reads


def run_isoseq_cluster(ccs_reads):

    log.info("ccs_reads " + ccs_reads)
    log.info('lima ...')
    demux_reads = "%s.demux.bam" % prefix
    cmd = 'lima --isoseq --dump-clips -j %s %s %s %s' % (n_threads,ccs_reads,primers,demux_reads)
    #log('lima',"cmd=%s" % cmd)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    #err_test(std_out, std_err, exit_code, cmd)

    demux_reads_p = "%s.demux.primer_5p--primer_3p.bam" % prefix
    unpolished_f = "%s.unpolished.flnc.bam" % prefix
    cmd = 'isoseq3 refine --require-polya %s %s %s' % (demux_reads_p,primers,unpolished_f)
    std_out, std_err, exit_code = run_cmd(cmd, log)
#isoseq3 refine --require-polya CZYTW.demux.primer_5p--primer_3p.bam ~/src/iso-seq3/primers.fasta CZYTW.unpolished.flnc.bam

#isoseq3 cluster -j 16 CZYTW.unpolished.flnc.bam CZYTW.unpolished.bam --split-bam 12

    
    log.info('isoseq3 cluster ...')
    #demux_reads_p="%s.demux.primer_5p--primer_3p.bam" % prefix
    #unpolished="%s.unpolished.flnc.bam" % prefix
    unpolished = "%s.unpolished.bam" % prefix
    cmd = 'isoseq3 cluster %s %s --split-bam %s -j %s ' % (unpolished_f,unpolished,n_split,n_threads)
    std_out, std_err, exit_code = run_cmd(cmd, log)

def run_isoseq_polish():
    log.info('isoseq3 polish, generate array script...')
    f = open("isoseq.polish.sh","w")
    f.write('#!/bin/bash -l\nexport j=`printf %06d $SLURM_ARRAY_TASK_ID`\nexport JOBNUM=$j\numask 2\nset -ex\n');
    f.write('./isoseq.polish.$SLURM_ARRAY_TASK_ID $JOBNUM\n')
    f.close()
    for i in range(0,n_split):
        f = open("isoseq.polish.%s" % str(i+1),"w")
        polished = "%s.polished.%s.bam" % (prefix,str(i))
        unpolished = "%s.unpolished.%s.bam" % (prefix,str(i))
        cmd = 'isoseq3 polish -j %s %s %s %s\n' % (n_threads,unpolished,subreads,polished)
        f.write('#!/bin/bash -l\n'+'source activate isoseq3_env\n'+'set -ex\n'+cmd)
        f.close()
        os.chmod("isoseq.polish.%s" % str(i+1),0o755)
    #sys.exit(1)
    cmd = 'sbatch -c %s -t 4:00:00 --parsable --array=1-%s isoseq.polish.sh' % (n_threads,n_split)
    #cmd = 'sbatch -c %s --time=24:00:00 -A gtrnd --qos=genepool --parsable --array=1-%s isoseq.polish.sh' % (n_threads,n_split)
    std_out, std_err, exit_code = run_cmd(cmd, log)

    jobid=std_out.rstrip()
    log.info("jobid=%s" % (jobid))
    done=0
    while done == 0:
        cmd = "sacct -b -j %s > polish.sacct.out" % jobid
        std_out, std_err, exit_code = run_cmd(cmd, log)
        cmd1 = "sacct -b -j %s |awk 'BEGIN{n=0}{if(/COMPLETED/) n=n+1 }END{print n}'" % jobid
        cmd2 = "sacct -b -j %s |awk 'BEGIN{n=0}{if(/FAILED/) n=n+1 }END{print n}'" % jobid
        std_out1, std_err, exit_code = run_cmd(cmd1,log)
        std_out2, std_err, exit_code = run_cmd(cmd2,log)
        std_out1 = std_out1.strip()
        std_out2 = std_out2.strip()
        print std_out1, std_out2

        if int(std_out1) == 2*int(n_split):
            done = 1
        elif int(std_out2) >0:
            log.error("polish failed, exit...")
            cmd = "scancel %s" % jobid
            std_out, std_err, exit_code = run_cmd(cmd, log)
            sys.exit(1)
        else:
            log.info("sleep 15*60 seconds")
            time.sleep(5*60)

    cmd = "gzip -f -d *.gz"
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "cat {0}.polished.*.hq.fastq > {0}.hq.fastq".format(prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "cat {0}.polished.*.hq.fasta > {0}.hq.fasta".format(prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "cat {0}.polished.*.lq.fastq > {0}.lq.fastq".format(prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "cat {0}.polished.*.lq.fasta > {0}.lq.fasta".format(prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "ls -1d {0}.polished.*.bam > list".format(prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "bamtools merge -list list -out {0}.polished.bam".format(prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "rm -f {0}.polished.*.bam  {0}.polished.*.bam.pbi {0}.polished.*.*q.fastq {0}.polished.*.*q.fasta".format(prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)


def run_minimap2(reference):
    log.info('run minimap2')
    fastq = '%s.hq.fastq' % prefix
    sam_sort = "%s.hq_isoforms.fastq.sorted.sam" % prefix
    cmd = 'shifter --image=robegan21/minimap2:2.10 minimap2 -t %s -ax splice -uf --secondary=no -C5 %s %s | sort -k 3,3 -k 4,4n > %s ' % (n_threads, reference, fastq, sam_sort)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    log.info("minimap2 complete")
    
def post_isoseq():
    log.info("post_isoseq ...")
    #cmd='%s/isoseq_post_process.sh %s %s' % (my_dir,prefix,primers)
    cmd = "isoseq3_make_classify_report.py  --flnc_bam {0}.unpolished.flnc.bam {0}.demux.lima.clips {1}".format(prefix, primers)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "isoseq3_make_cluster_report.py {0}.polished.bam".format(prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "collapse_isoforms_by_sam.py  --input {0}.hq.fastq --fq -s {0}.hq_isoforms.fastq.sorted.sam --dun-merge-5-shorter -o {0}.isoform".format(prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "get_abundance_post_collapse.py {0}.isoform.collapsed cluster_report.csv".format(prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "gffread {0}.isoform.collapsed.gff -T -o {0}.isoform.collapsed.gtf".format(prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "filter_away_subset.py {0}.isoform.collapsed".format(prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "gffread {0}.isoform.collapsed.filtered.gff -T -o {0}.isoform.collapsed.filtered.gtf".format(prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    log.info("post_isoseq complete")

  
def run_sqanti(reference,ann):
    log.info("sqanti ...")
    if not os.path.isdir("sqanti_out"):
         os.makedirs("sqanti_out")
    #cmd = "shifter --image=mjblow/sqanti:v1 sqanti_qc.py -d sqanti_out -g {0}.isoform.collapsed.gtf {1} {2} > sqanti.log 2>&1 ".format(prefix, ann, reference)
    cmd = "shifter --image=mjblow/sqanti:v1 sqanti_qc.py -d sqanti_out -g {0}.isoform.collapsed.filtered.gtf {1} {2} > sqanti.log 2>&1 ".format(prefix, ann, reference)
    std_out, std_err, exit_code =  run_cmd(cmd, log)
    log.info("sqanti complete")

def gather_stats():
    log.info("gather stats ...")
    dat = pd.DataFrame(index=[prefix])
  #  '''
    cmd = '{0} flagstat {1} > {1}.flagstat'.format(samtools,subreads)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = '{0} flagstat {1} > {1}.flagstat'.format(samtools, ccs_reads)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "%s view %s | awk -F'\t' '{print length($10)}' > %s.length.txt" % (samtools, subreads, subreads)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "%s view %s | awk -F'\t' '{print length($10)}' > %s.length.txt" % (samtools, ccs_reads, ccs_reads)
    std_out, std_err, exit_code = run_cmd(cmd, log)
   # '''
    bindir = os.path.dirname(os.path.realpath(__file__))
    cmd = 'shifter --volume=%s:/home/docker --image=docker:rmonti/r-essentials %s %s %s' % (bindir, prefix, subreads, ccs_reads)

    for line in open('{0}.flagstat'.format(subreads), 'r'):
        if 'in total' in line:
            dat['subreads_num'] = int(line.split()[0])

    for line in open('{0}.flagstat'.format(ccs_reads), 'r'):
        if 'in total' in line:
            dat['ccs_num'] = int(line.split()[0])
    
    cmd = '{0} view {1}.unpolished.flnc.bam |wc -l'.format(samtools, prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    dat['flnc_num'] = int(std_out.strip())

    cmd = 'grep -c ">" {0}.hq.fasta'.format(prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    dat['hq_isoform_num'] = int(std_out.strip())

    # parse mapping results
    if os.path.isfile("{0}.hq_isoforms.fastq.sorted.sam".format(prefix)):
        cmd = "awk -F'\t' '{if($2==0 || $2==16) print $1}' %s.hq_isoforms.fastq.sorted.sam |sort|uniq |wc -l" % (prefix)
        std_out, std_err, exit_code = run_cmd(cmd, log)
        dat['mapping_rate'] = int(std_out.strip())/float(dat['hq_isoform_num'])
        # get the isoform number based on mapping results    
        cmd = 'grep -wc "transcript" %s.isoform.collapsed.filtered.gff' % (prefix)
        std_out, std_err, exit_code = run_cmd(cmd, log)
        dat['mapped_isoform_num'] = int(std_out.strip())

    # parse sqanti output
    if gtf:
        cmd = "cut -f 9 %s" % gtf
        std_out, std_err, exit_code = run_cmd(cmd, log)
        tlist=[]
        for line in std_out.split('\n'):
            mobj = re.search(r'transcript_id "(.*?)";',line)
            if mobj:
                tlist.append(mobj.group(1))
        dat['ref_isoform_num'] = len(list(set(tlist)))
        
        metric_list = ['antisense', 'full-splice_match', 'genic', 'genic_intron', 'intergenic', 'novel_in_catalog', 'novel_not_in_catalog']
        for metric in metric_list:
            dat[metric] = -1
        if os.path.isfile("sqanti_out/{0}.isoform.collapsed.filtered_classification.txt".format(prefix)):
            cmd='cut -f6 sqanti_out/%s.isoform.collapsed.filtered_classification.txt|sort|uniq -cd' % prefix
            std_out, std_err, exit_code = run_cmd(cmd, log)
            for line in std_out.split('\n'):
                if line:
                    mnum, metric = line.strip().split()[0:2]
                    if metric in metric_list:
                        dat[metric] = mnum

    dat.to_csv(path_or_buf=prefix+'_stats.txt',sep='\t')
    
#---------------------------------------------------------------------------------#
# Main Program

script_name = __file__

if __name__ == "__main__":
    
    # Prog details
    my_name = script_name 
    version = "1.0"
    UMASK = 0o022
    #my_dir=os.path.dirname(os.path.realpath(__file__) )
    
    # Parse options
    usage = "isoseq.py <subreads/ccs.bam> <primers.fasta> <prefix>\n"
    usage += "* %s, version %s\n" % (my_name, version)
    
    #Args
    parser = argparse.ArgumentParser(usage = usage,description='Program to run iso-seq3 pipeline and do annotaion comparison')
    parser.add_argument("subreads",help="subread bam data")
    parser.add_argument("primers",help="primers fasta file")
    parser.add_argument("prefix",help="prefix name")
    parser.add_argument("-r","--reference",help="reference genome fasta file for post-processing",required=False)
    parser.add_argument("-c","--ccs",help="ccs reads", required=False)
    parser.add_argument("-g","--gtf",help="reference annotaion gtf file for running sqanti",required=False)
    parser.add_argument("-o", "--out_dir", dest="out_dir", help="output directory. default=./", default=".")
    parser.add_argument("-n", "--num_threads",help="number of threads used; default=16",default=16,type=int )
    parser.add_argument("-pl", "--print-log", dest="print_log", default=False, action="store_true", help="print log to screen")
    args=parser.parse_args()
    

    #setup
    setup()
        
    #run ccs on subreads
    if not args.ccs:
        log.info("CCS reads do not exit, start to generate it")
        ccs_reads=run_CCS()

#    '''
    # run iso_seq3 pipeline
    run_isoseq_cluster(ccs_reads)
    run_isoseq_polish()
 #   '''
                                     
    ###post isoseq_analysis
    if args.reference:
        log.info("reference files exits, strat to run minimap2 and squanti")
        run_minimap2(reference)
        post_isoseq()
        if args.gtf:
            run_sqanti(reference,gtf)
    #'''

    gather_stats()
            
    #complete, exit nicely
    log.info('iso-seq3 pipeine complete...')
    sys.exit(0)
