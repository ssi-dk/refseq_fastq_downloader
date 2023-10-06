import os
import pandas as pd
import argparse
import re
import subprocess
parser = argparse.ArgumentParser()
parser.add_argument("--input_file", "-i", default = 'gcf_accessions.txt', help = "text file with gcf accessions")
parser.add_argument('--output_dir', '-o', default = 'reads', help = 'directory to store the subdirectories')
parser.add_argument('--chunk_size', '-chunk', default = 1000, help = 'chunk_size', type=int)
args = parser.parse_args()
gcfs = []
with open(args.input_file) as input:
    lines = input.readlines()
    for i in lines:
        gcfs.append(i.strip())
if not os.path.isdir(args.output_dir):
    os.makedirs(args.output_dir)

chunk_size = args.chunk_size
chunks = list(range(0, len(gcfs)+chunk_size, chunk_size))

# for testing
#chunks = chunks[0]
for i in range(len(chunks)-1):
    idxs = [chunks[i], min(chunks[i+1], len(gcfs))]
    idxs_str = ','.join([str(i) for i in idxs])
    cmd = f'python run_fastq_dl_on_subset.py -i {args.input_file} -o {args.output_dir} -j {idxs_str}'
    slurm_cmd = f"sbatch -D . -c 4 --mem=10G -J getting_refseq_fasta_{idxs[0]}_{idxs[1]} -p fat02 --wrap=\'{cmd}\'"
    print(slurm_cmd)
    process = subprocess.Popen(slurm_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, env=os.environ, encoding='utf-8')
    process_out, process_err = process.communicate()
    print(process_out, process_err)

