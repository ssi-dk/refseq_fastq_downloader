import os
import pandas as pd
import argparse
import re
import subprocess
parser = argparse.ArgumentParser()
parser.add_argument("--input_file", "-i", default = 'gcf_accessions.txt', help = "text file with gcf accessions")
parser.add_argument('--output_dir', '-o', default = 'reads', help = 'directory to store the subdirectories')
parser.add_argument("--accession_idxs", "-j", default = None, help = 'indexes to slice accessions')
args = parser.parse_args()
gcfs = []
with open(args.input_file) as input:
    lines = input.readlines()
    for i in lines:
        gcfs.append(i.strip())

if args.accession_idxs != None:
    idxs = [int(i) for i in args.accession_idxs.split(',')]
    gcfs = gcfs[idxs[0]:idxs[1]]
print('running on ', len(gcfs), 'refseq accessions')
out = args.output_dir
for i in gcfs:
    cmd = f'python /srv/data/MPV/LEBC/dl_refseq_fastq/download_refseq_fastq.py -refseq {i} -o {out} --download_fastq'
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, env=os.environ, encoding='utf-8')
    process_out, process_err = process.communicate()
    print(process_out, process_err)
