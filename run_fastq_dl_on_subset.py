import os
import pandas as pd
import argparse
import re
import subprocess
import requests
import json
import xml.etree.ElementTree as ET
import sys
import hashlib
parser = argparse.ArgumentParser()
parser.add_argument("--input_file", "-i", help = "plain text file with newline separated gcf accessions OR tsv with 'GCF' column")
parser.add_argument('--output_dir', '-o', default = 'reads', help = 'directory to store the subdirectories')
parser.add_argument("--download_fastq", "-d", help="download fastq file", action="store_true", default=False)
parser.add_argument("--accession_idxs", "-j", default = None, help = 'indexes to slice accessions')
args = parser.parse_args()

def download_refseq_fastq(refseq_accession, output_folder, download_fastq):
    dataset_report = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{refseq_accession}/dataset_report"
    report_request = requests.get(dataset_report)
    status = report_request.status_code
    refseq_folder = os.path.join(output_folder, refseq_accession)
    if status == 200:
        report_json = report_request.content.decode()
        report_dict = json.loads(report_json)
        ena_accession = report_dict['reports'][0]['assembly_info']['biosample']['accession']
        ena_xml = f"https://www.ebi.ac.uk/ena/browser/api/xml/{ena_accession}?includeLinks=false"
        ena_request = requests.get(ena_xml)
        #print(ena_request.status_code)
        if int(ena_request.status_code) == 404:
            print(f'{refseq_accession},{ena_accession},no xml found')
            return
        #print(ena_request.content.decode())
        root = ET.fromstring(ena_request.content.decode())
        sample = root.find('SAMPLE')
        sample_dict = sample.attrib
        links = sample.find('SAMPLE_LINKS')
        sample_link = links.find('SAMPLE_LINK')
        xref_link = sample_link.find('XREF_LINK')
        link_id = xref_link.find('ID')
        fastq_locations = link_id.text
        ena_fastq_request = requests.get(fastq_locations)
        ena_fastq_bytes = len(ena_fastq_request.content)
        if ena_fastq_bytes > 46: # size of the dataframe header
            if not os.path.isdir(refseq_folder):
                os.mkdir(refseq_folder)
            with open(os.path.join(refseq_folder, 'ena_fastq.tsv'), 'w') as output:
                output.write(ena_fastq_request.content.decode())
            if download_fastq:
                ena_link_df = pd.read_csv(os.path.join(refseq_folder, 'ena_fastq.tsv'), sep='\t')
                ena_link_df = ena_link_df.dropna() # missing values will crash the code below
                #print(ena_link_df)
                link_pairs = ena_link_df['fastq_ftp'].tolist()
                link_pairs = [i.split(';') for i in link_pairs]
                md5_pairs = ena_link_df['fastq_md5'].tolist()
                md5_pairs = [i.split(';') for i in md5_pairs]
                for i,j in zip(link_pairs, md5_pairs):
                    for k,m in zip(i,j):
                        cmd = f'wget \"{k}\" -P {refseq_folder}'
                        fastq_basename = os.path.basename(k)
                        print(f'downloading {refseq_accession},{ena_accession},{fastq_basename}')
                        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, env=os.environ, encoding='utf-8')
                        process_out, process_err = process.communicate()
                        with open(os.path.join(refseq_folder, fastq_basename),"rb") as f:
                            bytes = f.read() # read file as bytes
                            readable_hash = hashlib.md5(bytes).hexdigest()
                        if readable_hash != m:
                                print(f'{refseq_accession},{ena_accession},md5sums {m}, {readable_hash} dont match')
                        #print(process_out, process_err)
                    
        else:
            print(f'{refseq_accession},{ena_accession},no_reads')
    else:
        print(f'failed to get json for {refseq_accession}, {status}')


def download_multiple_refseq_fastqs(input_file, output_folder, download_fastq, accession_idxs):
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)
    if bool(re.search('tsv', input_file)):
        gcf_tsv = pd.read_csv(input_file, sep='\t')
        gcfs = gcf_tsv['GCF'].tolist()
    else:
        gcfs = []
        with open(input_file) as input:
            lines = input.readlines()
        for i in lines:
            gcfs.append(i.strip())

    if accession_idxs != None:
        idxs = [int(i) for i in accession_idxs.split(',')]
        gcfs = gcfs[idxs[0]:idxs[1]]
    print('running on ', len(gcfs), 'refseq accessions')
    for i in gcfs:
        download_refseq_fastq(i, output_folder, download_fastq)

def main(args):
    download_multiple_refseq_fastqs(args.input_file, args.output_dir, args.download_fastq, args.accession_idxs)

if __name__ == '__main__':
    main(args)
