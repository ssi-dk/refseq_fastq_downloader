import pandas as pd
import argparse
import requests
import json
import xml.etree.ElementTree as ET
import os
import subprocess
import sys
import hashlib
parser = argparse.ArgumentParser()
parser.add_argument("--refseq_accession", "-refseq", help = "gcf bla bla bla")
parser.add_argument("--output_location", "-o", help = "place to put a folder named after the accession with output files", default='.')
parser.add_argument("--download_fastq", "-d", help="download fastq file", action="store_true", default=False)
args = parser.parse_args()
refseq_accession = args.refseq_accession
dataset_report = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{refseq_accession}/dataset_report"
report_request = requests.get(dataset_report)
status = report_request.status_code
refseq_folder = os.path.join(args.output_location, refseq_accession)
if status == 200:
    report_json = report_request.content.decode()
    report_dict = json.loads(report_json)
    ena_accession = report_dict['reports'][0]['assembly_info']['biosample']['accession']
    ena_xml = f"https://www.ebi.ac.uk/ena/browser/api/xml/{ena_accession}?includeLinks=false"
    ena_request = requests.get(ena_xml)
    #print(ena_request.status_code)
    if int(ena_request.status_code) == 404:
        print(f'{refseq_accession},{ena_accession},no xml found')
        sys.exit()
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
        if args.download_fastq:
            ena_link_df = pd.read_csv(os.path.join(refseq_folder, 'ena_fastq.tsv'), sep='\t')
            ena_link_df = ena_link_df.dropna()
            print(ena_link_df)
            link_pairs = ena_link_df['fastq_ftp'].tolist()
            print(link_pairs)
            link_pairs = [i.split(';') for i in link_pairs]
            md5_pairs = ena_link_df['fastq_md5'].tolist()
            md5_pairs = [i.split(';') for i in md5_pairs]
            for i,j in zip(link_pairs, md5_pairs):
                if type(link_pairs) == float: # weird issue if there's a NaN value in the pandas dataframe, should probably just iterate over the pandas indexes in the future
                    continue
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
