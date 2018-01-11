#!/usr/bin/env python3

'''
Usage: check_pattern.py [options] -j JUNC (-b BAM | -g BIGWIG) <rs_pattern>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -j JUNC --junc=JUNC            RS junction file.
    -b BAM --bam=BAM               Alignment Bam file.
    -g BIGWIG --bigwig=BIGWIG      Alignment BigWig file.
    -p THREAD --thread=THREAD      Threads. [default: 5]
    --chrom-size=CHROM_SIZE        Path of chrom size file (required for -b).
    --bin-size=SIZE                Split bin size. [default: 1000]
    --bin-num=NUM                  Split bin mumber. [default: 10]
    --fold-change=fold             Fold change cutoff. [default: 2]
    --iteration-time=TIME          Iteration time. [default: 10000]
    --cutoff=CUTOFF                Cutoff of p-value. [default: 0.01]
    --remote                       Convert remote bam file.
    --stranded                     Stranded sequencing.
'''

import sys
import os.path
from collections import defaultdict
from multiprocessing import Pool
import pyBigWig
import numpy as np
from seqlib.ngs import bam_to_bedgraph
from seqlib.helper import run_command
from seqlib.path import which

__author__ = 'Xiao-Ou Zhang <xiaoou.zhang@umassmed.edu>'
__version__ = '0.0.1'


def check_pattern(options):
    # parse options
    strand_flag = options['--stranded']
    remote_flag = options['--remote']
    thread = int(options['--thread'])
    # check bigwig
    if options['--bam']:
        if which('bedGraphToBigWig') is not None:
            if strand_flag:
                plus_bg, minus_bg = bam_to_bedgraph(options['--bam'],
                                                    url=remote_flag,
                                                    stranded=True)
                plus_bw = bg_to_bw(plus_bg, options['--chrom-size'])
                minus_bw = bg_to_bw(minus_bg, options['--chrom-size'])
            else:
                bg = bam_to_bedgraph(options['--bam'], url=remote_flag)
                bw = bg_to_bw(bg, options['--chrom-size'])
        else:
            sys.exit('Could not find bedGraphToBigWig!')
    else:
        if os.path.isfile(options['--bigwig']):
            bw = options['--bigwig']
    # parse junction file
    if os.path.isfile(options['--junc']):
        junc_path = options['--junc']
    else:
        sys.exit('Wrong junc file: %s' % options['--junc'])
    # check pattern
    junc_info = defaultdict(list)
    p = Pool(thread)
    result = []
    with open(junc_path, 'r') as junc_f:
        for junc in junc_f:
            info = junc.split()
            chrom = info[1]
            strand = info[4]
            rs_site = info[6].split('|')[0]
            rs_info = '\t'.join([chrom, rs_site, strand])
            if rs_info not in junc_info:
                if strand_flag:
                    if strand == '+':
                        result.append(p.apply_async(cal_ratio,
                                                    args=(plus_bw, rs_info,
                                                          options,)))
                    else:
                        result.append(p.apply_async(cal_ratio,
                                                    args=(minus_bw, rs_info,
                                                          options,)))
                else:
                    result.append(p.apply_async(cal_ratio, args=(bw, rs_info,
                                                                 options,)))
            junc_info[rs_info].append(info)
        p.close()
        p.join()
    with open(options['<rs_pattern>'], 'w') as out:
        for r in result:
            pvalue, fold, rs_info = r.get()
            if pvalue is not None:
                for junc in junc_info[rs_info]:
                    out.write('\t'. join(junc))
                    out.write('\t%f\t%f\n' % (fold, pvalue))


def bg_to_bw(bg, chrom_size):
    if not os.path.isfile(chrom_size):
        sys.exit('No chrom size file: %s!' % chrom_size)
    prefix = os.path.splitext(os.path.split(bg)[-1])[0]
    bw = prefix + '.bw'
    command = 'bedGraphToBigWig %s %s %s' % (bg, chrom_size, bw)
    run_command(command, 'Could not convert bg to bw!')
    return bw


def cal_ratio(bw, rs_info, options):
    # parse options
    bin_size = int(options['--bin-size'])
    bin_num = int(options['--bin-num'])
    fold_change = float(options['--fold-change'])
    iteration_time = int(options['--iteration-time'])
    cutoff = float(options['--cutoff'])
    # parse info
    bwf = pyBigWig.open(bw)
    chrom, site, strand = rs_info.split()
    site = int(site)
    # split bins
    region1, region2 = [], []
    for i in range(-bin_num, bin_num):
        if i < 0:
            signal = fetch_signal(bwf, chrom, site, bin_size, i)
            region1.append(signal)
        else:
            signal = fetch_signal(bwf, chrom, site, bin_size, i)
            region2.append(signal)
    region1 = np.array(region1)
    region2 = np.array(region2)
    # calculate difference
    d = (np.median(region1) + 0.1) / (np.median(region2) + 0.1)
    if strand == '+':
        peak_region = region2[0]  # potential as exon peak region
        if peak_region > np.mean(region2) + 3 * np.std(region2):  # as exon
            return None, None, rs_info
        if d <= 1 / fold_change:  # enough fold
            return permutation(region1, region2, d, 'less', iteration_time,
                               cutoff), d, rs_info
        else:  # not enough fold
            return None, None, rs_info
    else:
        peak_region = region1[-1]  # potential as exon peak region
        if peak_region > np.mean(region1) + 3 * np.std(region1):  # as exon
            return None, None, rs_info
        if d >= fold_change:  # enough fold
            return permutation(region1, region2, d, 'more', iteration_time,
                               cutoff), d, rs_info
        else:  # not enough fold
            return None, None, rs_info


def fetch_signal(bwf, chrom, site, bin_size, i):
    pos1 = site + i * bin_size
    pos2 = pos1 + bin_size
    signal = bwf.stats(chrom, pos1, pos2)[0]
    if signal is None:
        signal = 0
    return signal


def permutation(x, y, d, flag, ita, cutoff):
    pooled = np.hstack([x, y])
    dist = np.array(list(map(lambda i: rand(pooled, x.size, y.size),
                             range(ita))))
    if flag == 'more':
        count = len(np.where(dist >= d)[0])
    else:
        count = len(np.where(dist <= d)[0])
    pvalue = count * 1.0 / ita
    if pvalue < cutoff:
        return pvalue
    else:
        return None


def rand(pooled, xsize, ysize):
    np.random.shuffle(pooled)
    x = pooled[:xsize]
    y = pooled[-ysize:]
    return (np.median(x) + 0.1) / (np.median(y) + 0.1)


if __name__ == '__main__':
    from docopt import docopt
    check_pattern(docopt(__doc__, version=__version__))
