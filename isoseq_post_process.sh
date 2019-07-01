#!/bin/bash -l

# 1: prefix
# 2: primers file


function log {
    echo "$(date '+%C%y-%m-%d %T,%3N')    $1"
}

# exit on error function:
function eoe {
    if [ "$1" -ne 0 ]
        then
        echo "exit on error code $1"
        exit 1
    fi
}


log "start post processing"
#log "cmd=isoseq3_make_classify_report.py  --flnc_bam $1.unpolished.flnc.bam $1.demux.lima.clips $2"
#isoseq3_make_classify_report.py --flnc_bam $1.unpolished.flnc.bam $1.demux.lima.clips $2
#eoe "$?"

log "cmd=isoseq3_make_cluster_report.py $1.polished.bam"
isoseq3_make_cluster_report.py $1.polished.bam
eoe "$?"

log "cmd=collapse_isoforms_by_sam.py  --input $1.hq.fastq --fq -s $1.hq_isoforms.fastq.sorted.sam --dun-merge-5-shorter -o $1.isoform"
collapse_isoforms_by_sam.py --input $1.hq.fastq --fq -s $1.hq_isoforms.fastq.sorted.sam --dun-merge-5-shorter -o $1.isoform
eoe "$?"

log "cmd=get_abundance_post_collapse.py $1.isoform.collapsed cluster_report.csv"
get_abundance_post_collapse.py $1.isoform.collapsed cluster_report.csv
eoe "$?"

log "cmd=gffread $1.isoform.collapsed.gff -T -o $1.isoform.collapsed.gtf"
gffread $1.isoform.collapsed.gff -T -o $1.isoform.collapsed.gtf
eoe "$?"

log "complete post processing"
