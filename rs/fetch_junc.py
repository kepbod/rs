#!/usr/bin/env python3

'''
Usage: fetch_junc.py [options] -s SITE (-b BAM | -j JUNC) <rs_junc>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -s SITE --site=SITE            RS sites.
    -b BAM --bam=BAM               Alignment BAM file.
    -j JUNC --junc=JUNC            Junction BED file.
    --min-reads=READS              Minimum junction reads. [default: 2]
    --remote                       Fetch junction file from remote bam file.
    --stranded                     Stranded sequencing.
'''

import sys
import os.path
from seqlib.ngs import fetch_juncfile

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu'
__version__ = '0.0.1'


def fetch_junc(options):
    # parse options
    min_junc_read = int(options['--min-reads'])
    stranded_flag = options['--stranded']
    # parse junction file
    if options['--bam']:  # fetch junction from bam file
        junc_f = fetch_juncfile(options['--bam'], url=options['--remote'],
                                stranded=stranded_flag)
    else:  # direct use junction file
        if os.path.isfile(options['--junc']):
            junc_f = options['--junc']
        else:
            sys.exit('Wrong junction file: %s' % options['--junc'])
    junc = {}
    with open(junc_f, 'r') as f:
        for line in f:
            chrom, start, end, name, _, strand = line.split()[:6]
            s1 = int(start) + 10
            s2 = int(end) - 10
            if stranded_flag:
                junc_pos = '%s\t%d\t%d\t%s' % (chrom, s1, s2, strand)
            else:
                junc_pos = '%s\t%d\t%d\t+' % (chrom, s1, s2)
            junc_read = int(name.split('/')[-1])
            if junc_read >= min_junc_read:  # enough junction reads
                junc[junc_pos] = junc_read
    # parse rs site file
    if os.path.isfile(options['--site']):
        site_path = options['--site']
    else:
        sys.exit('Wrong site file: %s' % options['--site'])
    with open(site_path, 'r') as f, open(options['<rs_junc>'], 'w') as outf:
        for line in f:
            gene, chrom, start, end, strand, length, rs = line.split()
            if strand == '+':
                rs_junc = chrom + '\t' + start + '\t%s'
            else:
                rs_junc = chrom + '\t%s\t' + end
            if stranded_flag:
                rs_junc += '\t%s' % strand
            else:
                rs_junc += '\t+'
            for rs_info in rs.split(','):
                rs_site = rs_info.split('|')[0]
                rs_junc_pos = rs_junc % rs_site
                if rs_junc_pos in junc:
                    outf.write('\t'.join([gene, chrom, start, end, strand,
                                          length, rs_info,
                                          str(junc[rs_junc_pos])]))
                    outf.write('\n')


if __name__ == '__main__':
    from docopt import docopt
    fetch_junc(docopt(__doc__, version=__version__))
