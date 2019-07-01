# New pipeline for PacBio Sequel II data
# Based off of JGI Jie Guo's isoseq.py pipeline,
#   and uses much of her code for gather_stats().
# The functions in PostProcessing.generate_reports()
#   were written by Elizabeth Tseng at PacBio. 

import os
import sys
import argparse
import re
import subprocess

import numpy as np
import pandas as pd

from jgi_common_funcs import run_cmd, get_logger

class CCS:

    def __init__(self, subreads, prefix, out_dir, split_num=0):
        if not subreads.endswith('.bam'):
            log.error("Subreads file must be in bam format!")
            sys.exit(1)
        self.subreads = subreads
        self.prefix = prefix
        self.out_dir = out_dir
        if split_num == 0:
            size = os.path.getsize(subreads)
            log.info("Size of subreads: %dGB" % (size/10**9))
            if size <= 40*(10**9):
                self.split_num = 1
            else:
                chunks = size/10**9/40
                log.info("Will break ccs into %d jobs" % (chunks))
                self.split_num = chunks
        else:
            self.split_num = split_num

    def get_zmw_ranges(self):
        zmw_file = '%s_ZMWs.txt' % self.prefix
        bamsieve_cmd = 'bamSieve --show-zmws %s > %s' % (self.subreads, zmw_file)
        run_cmd(bamsieve_cmd)
        with open(zmw_file) as zmws:
            total_lines=0
            for line in zmws:
                total_lines += 1
        region_size = total_lines / self.split_num
        with open(zmw_file) as zmws:
            range_start = int(zmws.readline().strip())
            count=1
            ranges = []
            range_end = None
            for line in zmws:
                count+=1
                if range_end:
                    ranges.append((range_start, range_end))
                    range_start = int(line.strip())
                    range_end = None
                    continue
                if count % region_size == 0:
                    range_end = int(line.strip())
                if count == total_lines:
                    range_end = int(line.strip())
                    ranges.append((range_start, range_end))
        # The very last range is a very small region of leftover zmws.
        # I want to include them in the second to last range,
        # So I remove the last two ranges, and add a new last range
        # that uses the start of the second to last range and the end of the last range
        last_start, last_end = ranges.pop()
        scd_last_start, scd_last_end = ranges.pop()
        ranges.append((scd_last_start, last_end))
        # Writing regions as a data file to input into ccs array cmd
        array_data_text = ""
        region_num = 0
        for region in ranges:
            region_num += 1
            array_data_text += '%d-%d\t%d\n' % (region[0], region[1], region_num)
        write_file('%s_ccs_array_data.txt' % self.prefix, filetype='text', text=array_data_text)

    def write_ccs_command(self):
        if self.split_num >= 2:
            wrapper_code = "ZMWS=$1\nREGION_NUM=$2\n\n"
            wrapper_code += 'date\n'
            wrapper_code += "ccs --minReadScore=0.65 --minZScore=-3.5 --numThreads=32 --minSnr=4 "
            wrapper_code += "--maxLength=21000 --minLength=10 --minPasses=3 --minPredictedAccuracy=0.98 "
            wrapper_code += '--maxDropFraction=0.34 --reportFile=ccs_report_"${REGION_NUM}".txt '
            wrapper_code += '--zmws="${ZMWS}" %s "%s.ccs_${REGION_NUM}.bam"\n' % (self.subreads, self.prefix)
            wrapper_code += 'date\n'
            write_file('%s_ccs_array_wrapper.sh' % self.prefix, filetype='shell', text=wrapper_code)
            self.ccs_file = '%s_ccs_array_cmd.sl' % self.prefix
            write_file(self.ccs_file, filetype='array', array_data='%s_ccs_array_data.txt' % self.prefix,
                       array_wrapper='%s_ccs_array_wrapper.sh' % self.prefix, qos='genepool', tasks_per_node=32,
                       time='20:00:00', array_range='1-%d' % self.split_num,
                       slurm_out='{0}/slurm_output/ccs_%a.out'.format(self.out_dir))
        else:
            ccs_code = "date\n"
            ccs_code += "ccs --minReadScore=0.65 --minZScore=-3.5 --numThreads=32 --minSnr=4 "
            ccs_code += "--maxLength=21000 --minLength=10 --minPasses=3 --minPredictedAccuracy=0.98 "
            ccs_code += '--maxDropFraction=0.34 --reportFile=ccs_report.txt '
            ccs_code += '%s %s.ccs.merged.bam\n' % (self.subreads, self.prefix)
            ccs_code += 'date\n'
            self.ccs_file = '%s_ccs_cmd.sh' % self.prefix
            write_file(self.ccs_file, filetype='sbatch', text=ccs_code, qos='genepool',
                       tasks_per_node=1, time='20:00:00', slurm_out='%s/slurm_output/ccs.out' % self.out_dir)
        self.ccs_cmd = 'sbatch --parsable %s' % self.ccs_file


    def generate_ccs(self):
        if self.split_num >= 2:
            log.info("Calculating ZMW ranges...")
            self.get_zmw_ranges()
        log.info("Writing ccs command...")
        self.write_ccs_command()
        log.info("Submitting ccs commnd...")
        std_out, std_err, exit_code = run_cmd(self.ccs_cmd)
        ccs_jobid = std_out.rstrip()
        log.info("CCS cmd: %s" % self.ccs_cmd)
        log.info("CCS jobid: %s" % ccs_jobid)
        ccs_merge = False
        if self.split_num >= 2:
            ccs_merge = True
        return ccs_jobid, ccs_merge
        

class IsoSeq3:

    def __init__(self, subreads, primers, prefix, out_dir, ccs_jobid, ccs_merge):
        self.subreads = subreads
        self.primers = primers
        self.prefix = prefix
        self.out_dir = out_dir
        self.ccs_jobid = ccs_jobid
        self.ccs_merge = ccs_merge
        self.ccs = '%s.ccs.merged.bam' % prefix
        self.demux = "%s.demux.bam" % prefix
        self.demux_reads_p = '%s.demux.primer_5p--primer_3p.bam' % prefix
        self.unpolished_flnc = '%s.unpolished.flnc.bam' % prefix
        self.unpolished = '%s.unpolished.bam' % prefix

    def cluster(self):
        # First, merge together the multiple ccs files generated 
        #    from the ccs commands if it needed to be parallelized.
        # run Lima, the PacBio barcode demultiplexer,
        # and then the isoseq3 refine cmd to generate Full-length non-chimeric (FLNC) reads.
        # Then run isoseq3 polish command to get HQ isoforms.
        pipeline_prefix = self.prefix.split('/')[-1]
        cluster_code = 'date\n'
        if self.ccs_merge == True:
            cluster_code += 'python isoseq_pipeline.py --ccs_merge --out_dir %s --prefix %s\n' % (self.out_dir, pipeline_prefix)
        cluster_code += 'lima --isoseq --dump-clips -j 32 %s %s %s\ndate\n' % (self.ccs, self.primers, self.demux)
        cluster_code += 'isoseq3 refine --require-polya %s %s %s\ndate\n' % (self.demux_reads_p, self.primers,
                                                                             self.unpolished_flnc)
        cluster_code += 'isoseq3 cluster %s %s --split-bam 36 -j 32\ndate\n' % (self.unpolished_flnc,
                                                                                self.unpolished)
        self.cluster_file = '%s_cluster_cmd.sh' % self.prefix
        write_file(self.cluster_file, filetype='sbatch', text=cluster_code, qos='genepool',
                   tasks_per_node=32, time='2:00:00', slurm_out='%s/slurm_output/cluster.out' % self.out_dir)
        if self.ccs_jobid != 0:
            cluster_cmd = 'sbatch --parsable --dependency=afterok:%s %s' % (self.ccs_jobid, self.cluster_file)
        else:
            cluster_cmd = 'sbatch --parsable %s' % self.cluster_file
        std_out, std_err, exit_code = run_cmd(cluster_cmd)
        self.cluster_jobid = std_out.rstrip()
        log.info("Cluster cmd: %s" % cluster_cmd)                                                                                                

    def polish(self):
        wrapper_unpolished = '"%s.unpolished.${NUM}.bam"' % self.prefix
        wrapper_polished = '"%s.polished.${NUM}.bam"' % self.prefix
        wrapper_code = 'NUM=$(expr $1 - 1)\n\n'
        wrapper_code += 'date\n'
        wrapper_code += 'isoseq3 polish -j 32 %s %s %s' % (wrapper_unpolished, self.subreads,
                                                           wrapper_polished)
        wrapper_code += '\ndate\n'
        self.wrapper_file = '%s_polish_array_wrapper.sh' % self.prefix
        write_file(self.wrapper_file, filetype='shell', text=wrapper_code)
        # slightly different strategy for this array,
        # so writing the code explicitly and calling write_file() as a slurm file.
        array_code = '%s "$SLURM_ARRAY_TASK_ID"\n' % self.wrapper_file
        self.array_file = '%s_polish_array_cmd.sl' % self.prefix
        write_file(self.array_file, filetype='sbatch', text=array_code, qos='genepool',
                   tasks_per_node=32, time='2:00:00', array_range='1-36',
                   slurm_out='{0}/slurm_output/polish_%a.out'.format(self.out_dir))
        polish_cmd = 'sbatch --parsable --dependency=afterok:%s %s' % (self.cluster_jobid, self.array_file)
        std_out, std_err, exit_code = run_cmd(polish_cmd)
        self.polish_jobid = std_out.rstrip()
        log.info("Polish cmd: %s" % polish_cmd)

    def isoseq3_pipeline(self):
        log.info("Clustering ccs data...")
        self.cluster()
        log.info("Cluster jobid: %s" % self.cluster_jobid)
        log.info("Polishing unpolished bam files...")
        self.polish()
        log.info("Polish jobid: %s" % self.polish_jobid)
        return self.polish_jobid

class PostProcessing:

    def __init__(self, subreads, reference, gtf, primers, prefix, out_dir, polish_jobid, pipeline_dir):
        self.subreads = subreads
        self.reference = reference
        self.gtf = gtf
        self.primers = primers
        self.prefix = prefix
        self.out_dir = out_dir
        self.polish_jobid = polish_jobid
        self.pipeline_dir = pipeline_dir
        self.ccs = '%s.ccs.merged.bam' % prefix
        self.hq = '%s.hq.fastq' % prefix
        self.hq_sorted = '%s.hq_isoforms.fastq.sorted.sam' % prefix

    def mapping(self):
        # first: combine polished files into fasta/q files,
        #    and remove all of the polished and unpolished bam files
        # Then, Map with minimap2
        mapping_code = 'gzip -f -d %s/*.gz\n' % self.out_dir
        file_exts = ('hq.fastq', 'hq.fasta', 'lq.fastq', 'lq.fasta')
        for file_ext in file_exts:
            mapping_code += 'cat {0}.polished.*.{1} > {0}.{1}\n'.format(self.prefix, file_ext)
        polished_list = '%s.polished_list.txt' % self.prefix
        mapping_code += "ls -ld %s.polished.*.bam | cut -f 9 -d ' ' > %s\n" % (self.prefix, polished_list)
        mapping_code += 'bamtools merge -list %s -out %s.polished.bam\n' % (polished_list, self.prefix)
        mapping_code += 'rm -f {0}.polished.*.bam* {0}.polished.*.*q.fast*\n'.format(self.prefix)
        mapping_code += 'shifter --image=robegan21/minimap2:2.10 minimap2 -t 32 -ax splice -uf '
        mapping_code += '--secondary=no -C5 %s %s | sort -k 3,3 -k 4,4n > %s' % (self.reference, self.hq,
                                                                                 self.hq_sorted)
        self.mapping_file = '%s_mapping_cmd.sh' % self.prefix
        write_file(self.mapping_file, filetype='sbatch', text=mapping_code, qos='genepool',
                   tasks_per_node=32, time='04:00:00', slurm_out='%s/slurm_output/mapping.out' % self.out_dir)
        if self.polish_jobid != 0:
            mapping_cmd = 'sbatch --parsable --dependency=afterok:%s %s' % (self.polish_jobid, self.mapping_file)
        else:
            mapping_cmd = 'sbatch --parsable %s' % self.mapping_file
        std_out, std_err, exit_code = run_cmd(mapping_cmd)
        self.mapping_jobid = std_out.rstrip()
        log.info("Mapping cmd: %s" % mapping_cmd)

    def collapse(self):
        # Collapse with cupcake
        collapse_code = 'python %s/collapse_isoforms_by_sam.py --input %s.hq.fastq' % (self.pipeline_dir, self.prefix)
        collapse_code += '--fq -s {0}.hq_isoforms.fastq.sorted.sam --dun-merge-5-shorter -o {0}.isoform'.format(self.prefix)
        self.collapse_file = "%s_collapse_cmd.sh" % self.prefix
        write_file(self.collapse_file, filetype='sbatch', text=collapse_code, qos='genepool',
                   tasks_per_node=32, time='04:00:00', slurm_out='%s/slurm_output/collapse.out' % self.out_dir)
        collapse_cmd = 'sbatch --parsable --dependency=afterok:%s %s' % (self.mapping_jobid, self.collapse_file)
        std_out, std_err, exit_code = run_cmd(collapse_cmd)
        self.collapse_jobid = std_out.rstrip()
        log.info("Collapse cmd: %s" % collapse_cmd)

    def generate_reports(self):
        report_code = 'python %s/isoseq3_make_classify_report.py --flnc_bam ' % self.pipeline_dir
        report_code += '{0}.unpolished.flnc.bam {0}.demux.lima.clips {1}\n'.format(self.prefix, self.primers)
        report_code += 'python %s/isoseq3_make_cluster_report.py %s.polished.bam\n' % (self.pipeline_dir, self.prefix)
        report_code += 'python %s/collapse_isoforms_by_sam.py --input %s ' % (self.pipeline_dir, self.hq)
        report_code += ' --fq -s %s --dun-merge-5-shorter -o %s.isoform\n' % (self.hq_sorted, self.prefix)
        report_code += 'python %s/get_abundance_post_collapse.py ' % self.pipeline_dir
        report_code += '%s.isoform.collapsed cluster_report.csv\n' % self.prefix
        report_code += 'gffread {0}.isoform.collapsed.gff -T -o {0}.isoform.collapsed.gtf\n'.format(self.prefix)
        report_code += 'python %s/filter_away_subset.py %s.isoform.collapsed\n' % (self.pipeline_dir, self.prefix)
        report_code += 'gffread %s.isoform.collapsed.filtered.gff -T -o ' % self.prefix
        report_code += '%s.isoform.collapsed.filtered.gtf\n' % self.prefix
        self.report_file = '%s_report_cmd.sh' % self.prefix
        write_file(self.report_file, filetype='sbatch', text=report_code, qos='genepool',
                   tasks_per_node=32, time='04:00:00', slurm_out='%s/slurm_output/reports.out' % self.out_dir)
        report_cmd = 'sbatch --parsable --dependency=afterok:%s %s' % (self.collapse_jobid, self.report_file)
        std_out, std_err, exit_code = run_cmd(report_cmd)
        self.report_jobid = std_out.rstrip()
        log.info("Report cmd: %s" % report_cmd)

    def reference_comparison(self):
        # Compare to existing annotation w/ Sqanti
        sqanti_outdir = "%s/sqanti_out" % self.out_dir
        if not os.path.isdir(sqanti_outdir):
            os.mkdir(sqanti_outdir)
        prefix = self.prefix.split('/')[-1]
        sqanti_code = 'shifter --image=mjblow/sqanti:v1 sqanti_qc.py -d sqanti_out '
        sqanti_code += '-g %s.isoform.collapsed.filtered.gtf %s %s ' % (prefix, self.gtf, self.reference)
        sqanti_code += '> sqanti.log 2>&1 \n'
        self.sqanti_file = '%s_sqanti_cmd.sh' % self.prefix
        write_file(self.sqanti_file, filetype='sbatch', text=sqanti_code, qos='genepool',
                   tasks_per_node=32, time='04:00:00', slurm_out='%s/slurm_output/sqanti.out' % self.out_dir)
        sqanti_cmd = 'sbatch --parsable --dependency=afterok:%s %s' % (self.report_jobid, self.sqanti_file)
        std_out, std_err, exit_code = run_cmd(sqanti_cmd)
        self.sqanti_jobid = std_out.rstrip()
        log.info("Sqanti cmd: %s" % sqanti_cmd)

    def gather_stats(self):
        prefix = self.prefix.split('/')[-1]
        stats_code = 'python isoseq_pipeline.py --gather_stats --prefix %s ' % prefix
        stats_code += '--out_dir %s --gtf %s\n' % (self.out_dir, self.gtf)
        self.stats_file = '%s_gather_stats_cmd.sh' % self.prefix
        write_file(self.stats_file, filetype='sbatch', text=stats_code, qos='genepool',
                   tasks_per_node=32, time='06:00:00', slurm_out='%s/slurm_output/stats.out' % self.out_dir)
        stats_cmd = 'sbatch --parsable --dependency=afterok:%s %s' % (self.sqanti_jobid, self.stats_file)
        std_out, std_err, exit_code = run_cmd(stats_cmd)
        self.stats_jobid = std_out.rstrip()
        log.info("Stats cmd: %s" % stats_cmd)

    def postprocessing_pipeline(self):
        log.info("Mapping isoforms to reference genome...")
        self.mapping()
        log.info("Mapping jobid: %s" % self.mapping_jobid)
        log.info("Collapsing isoforms...")
        self.collapse()
        log.info("Collapsing jobid: %s" % self.collapse_jobid)
        log.info("Generating summary reports...")
        self.generate_reports()
        log.info("Report jobid: %s" % self.report_jobid)
        log.info("Comparing to existing annotation...")
        self.reference_comparison()
        log.info("Sqanti jobid: %s" % self.sqanti_jobid)
        log.info("Gathering stats...")
        self.gather_stats()
        log.info("Stats jobid: %s" % self.stats_jobid)

def write_file(filename, text=None, filetype='sbatch', qos='genepool_shared', 
               nodes=1, tasks_per_node=1, time='00:01:00', slurm_out=None,
               array_data=None, array_wrapper=None, array_range=None):
    if filetype not in ['text', 'shell', 'sbatch', 'array']:
        log.error("Cannot write filetype %s" % filetype)
        sys.exit(1)
    if filetype != 'array' and text == None:
        log.error("If filetype isn't array, must supply text for the file!")
        sys.exit(1)
    if filetype == 'array' and not all((array_range, os.path.exists(array_data), os.path.exists(array_wrapper))):
        log.error("If writing an array file, must supply array range plus data and wrapper files!")
        sys.exit(1)
    with open(filename, "w+") as f:
        if filetype != 'text':
            f.write("#!/bin/sh\n\n")
        if filetype in ['sbatch', 'array']:
            f.write("#SBATCH -A gtrnd\n#SBATCH -C Haswell\n")
            if array_range:
                f.write("#SBATCH --array=%s\n" % array_range)
            f.write("#SBATCH --qos=%s\n#SBATCH --nodes=%d\n" % (qos, nodes))
            f.write("#SBATCH --tasks-per-node=%d\n" % tasks_per_node)
            f.write("#SBATCH --time=%s\n\n" % time)
            if slurm_out:
                f.write("#SBATCH --output=%s\n" % slurm_out)
        if filetype != 'array':
            f.write(text)
        if filetype == 'array':
            f.write('ARGS=$(sed -n ${SLRUM_ARRAY_TASK_ID}p %s)\n\n' % array_data)
            f.write('%s ${ARGS}\n' % array_wrapper)
    run_cmd('chmod u+x %s' % filename)

def ccs_merge(prefix, out_dir):
    with open("%s.ccs_files.txt" % prefix, "w+") as outf:
        for root, subdirs, files in os.walk(out_dir):
            for file in files:
                if re.match('%s.ccs_[0-9]{1,2}.bam$' % prefix, file):
                    outf.write(file + '\n')
    log.info("Merging ccs bamfiles...")
    subprocess.check_output(['bamtools', 'merge', '-list', '%s.ccs_files.txt' % prefix,
                             '-out', '%s.ccs.merged.bam' % prefix])
    log.info("Indexing merged bamfile...")
    subprocess.check_output(['pbindex', '%s.ccs.merged.bam' % prefix])

def gather_stats(path_prefix, out_dir, subreads, gtf):
    log.info("gather stats ...")
    prefix = path_prefix.split('/')[-1]
    ccs_reads = '%s.ccs.merged.bam' % path_prefix
    dat = pd.DataFrame(index=[prefix])
    cmd = 'samtools flagstat {0} > {0}.flagstat'.format(subreads)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = 'samtools flagstat {0} > {0}.flagstat'.format(ccs_reads)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "%s view %s | awk -F'\t' '{print length($10)}' > %s.length.txt" % (samtools, subreads, subreads)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    cmd = "%s view %s | awk -F'\t' '{print length($10)}' > %s.length.txt" % (samtools, ccs_reads, ccs_reads)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    bindir = out_dir
    cmd = 'shifter --volume=%s:/home/docker --image=docker:rmonti/r-essentials %s %s %s' % (bindir, prefix, subreads, ccs_reads)

    for line in open('{0}.flagstat'.format(subreads), 'r'):
        if 'in total' in line:
            dat['subreads_num'] = int(line.split()[0])

    for line in open('{0}.flagstat'.format(ccs_reads), 'r'):
        if 'in total' in line:
            dat['ccs_num'] = int(line.split()[0])
    
    cmd = 'samtools view {0}.unpolished.flnc.bam | wc -l'.format(path_prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    dat['flnc_num'] = int(std_out.strip())

    cmd = 'grep -c ">" {0}.hq.fasta'.format(path_prefix)
    std_out, std_err, exit_code = run_cmd(cmd, log)
    dat['hq_isoform_num'] = int(std_out.strip())

    # parse mapping results
    if os.path.isfile("{0}.hq_isoforms.fastq.sorted.sam".format(path_prefix)):
        cmd = "awk -F'\t' '{if($2==0 || $2==16) print $1}' %s.hq_isoforms.fastq.sorted.sam |sort|uniq |wc -l" % (path_prefix)
        std_out, std_err, exit_code = run_cmd(cmd, log)
        dat['mapping_rate'] = int(std_out.strip())/float(dat['hq_isoform_num'])
        # get the isoform number based on mapping results    
        cmd = 'grep -wc "transcript" %s.isoform.collapsed.filtered.gff' % (path_prefix)
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
        if os.path.isfile("{0}/sqanti_out/{1}.isoform.collapsed.filtered_classification.txt".format(out_dir, prefix)):
            cmd='cut -f6 %s/sqanti_out/%s.isoform.collapsed.filtered_classification.txt|sort|uniq -cd' % (out_dir, prefix)
            std_out, std_err, exit_code = run_cmd(cmd, log)
            for line in std_out.split('\n'):
                if line:
                    mnum, metric = line.strip().split()[0:2]
                    if metric in metric_list:
                        dat[metric] = mnum

    dat.to_csv(path_or_buf=path_prefix+'_stats.txt',sep='\t')
    log.info("Analysis complete!")

def setup():
    correct_env = '$SCRATCH/Anaconda/IsoSeqEnv'
    std_out, std_err, exit_code = run_cmd('echo $CONDA_PREFIX')
    if not std_out.strip():
        print("No conda environment! Exiting...")
        sys.exit(1)
    elif 'IsoSeqEnv' not in std_out:
        print("Wrong conda environment! Exiting...")
        sys.exit(1)
    args = collect_args()
    global log
    log_name = "%s.log" % args.prefix
    if os.path.exists(log_name) and not any((args.ccs_merge, args.gather_stats)):
        os.remove(log_name)
    log = get_logger(__name__, log_name)
    return args

def collect_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--subreads', help="Path to subreads file")
    parser.add_argument('--gtf', help="Path to gtf file")
    parser.add_argument('--reference', help="Path to reference genome")
    parser.add_argument('--primers', help="Path to primers fasta file")
    parser.add_argument('--prefix', help="Prefix name")
    parser.add_argument('--out_dir', help="Output directory. Default ./", default='.')
    parser.add_argument('--resume_from', choices=['ccs', 'isoseq3', 'postprocessing'],
                         default='ccs',
                         help='Part of pipeline to start from, if some already completed')
    internal_use = parser.add_argument_group('internal_use')
    internal_use.add_argument('--ccs_merge', action='store_true',
                              help="Merges and index CCS files, then exits.")
    internal_use.add_argument('--gather_stats', action='store_true',
                              help="Gathers summary stats, then exits.")
    args = parser.parse_args()

    if args.ccs_merge:
        if not args.prefix:
            print("Prefix arg required with ccs_merge!")
            sys.exit(1)
    if args.gather_stats:
        if not all(args.prefix, args.subreads, args.gtf):
            print("Prefix, subreads, and gtf args required with ccs_merge!")
            sys.exit(1)
    if not all((args.subreads, args.gtf, args.reference, args.primers, args.prefix)):
        print("Subreads, gtf, reference, primers and prefix args must be given!")
        sys.exit(1)
    args.subreads, args.gtf, args.reference, args.primers, args.out_dir = (os.path.abspath(args.subreads),
                                                                           os.path.abspath(args.gtf),
                                                                           os.path.abspath(args.reference),
                                                                           os.path.abspath(args.primers),
                                                                           os.path.abspath(args.out_dir))
    for filepath in (args.subreads, args.gtf, args.reference, args.primers):
        if not os.path.exists(filepath):
            print("File %s does not exist! Exiting..." % filepath)
            sys.exit(1)
    if not os.path.isdir(args.out_dir):
        print("Dir %s does not exist! Exiting..." % args.out_dir)
        sys.exit(1)
    if not os.path.isdir("%s/slurm_output" % args.out_dir):
        os.mkdir("%s/slurm_output" % args.out_dir)
    args.prefix = '/'.join((args.out_dir, args.prefix))
    return args

if __name__ == '__main__':
    args = setup()
    pipeline_dir = os.path.dirname(__file__)
    log.info("Args: %s" % args)
    if args.ccs_merge:
        ccs_merge(args.prefix, args.out_dir)
        sys.exit(0)
    if args.gather_stats:
        gather_stats(args.prefix, args.out_dir, args.subreads, args.gtf)
        sys.exit(0)
    if args.resume_from == 'ccs':
        ccs_inst = CCS(args.subreads, args.prefix, args.out_dir)
        ccs_jobid, ccs_merge = ccs_inst.generate_ccs()
    else:
        log.info("Skipping ccs steps")

    if args.resume_from == 'ccs':
        isoseq3_inst = IsoSeq3(args.subreads, args.primers, args.prefix, args.out_dir, 
                               ccs_jobid, ccs_merge)
        polish_jobid = isoseq3_inst.isoseq3_pipeline()
    elif args.resume_from == 'isoseq3':
        log.info("Starting from isoseq3")
        isoseq3_inst = IsoSeq3(args.subreads, args.primers, args.prefix, args.out_dir, 
                               0, ccs_merge)
        polish_jobid = isoseq3_inst.isoseq3_pipeline()
    else:
        log.info("Skipping isoseq3 steps")
    if args.resume_from in ['ccs', 'isoseq3']:    
        postprocessing_inst = PostProcessing(args.subreads, args.reference, args.gtf, 
                                             args.primers, args.prefix, args.out_dir,
                                             polish_jobid, pipeline_dir)
    elif args.resume_from == 'postprocessing':
        log.info("Starting from postprocessing")
        postprocessing_inst = PostProcessing(args.subreads, args.reference, args.gtf, 
                                             args.primers, args.prefix, args.out_dir,
                                             0, pipeline_dir)

    postprocessing_inst.postprocessing_pipeline()
    log.info("All jobs submitted!")
    sys.exit(0)
