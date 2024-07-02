import os
import pandas as pd
import argparse
import re
import subprocess
import requests
import json
import sys
import hashlib
parser = argparse.ArgumentParser()
parser.add_argument("--input_file", "-i", help = "plain text file with newline separated gcf accessions OR tsv with 'GCF' column")
parser.add_argument('--output_tsv', '-o', default = 'asdf.tsv', help = 'output file with meta')
parser.add_argument("--accession_idxs", "-j", default = None, help = 'indexes to slice accessions')
args = parser.parse_args()

def download_refseq_meta(refseq_accession):
    dataset_report = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{refseq_accession}/dataset_report"
    report_request = requests.get(dataset_report)
    status = report_request.status_code
    if status == 200:
        report_content = report_request.content.decode()
        report_json = json.loads(report_content)
        if len(report_json.keys()) < 1:
            return None
        report_dict = report_json['reports'][0]
        assembly_info = report_dict.get('assembly_info')
        biosample = assembly_info.get('biosample')
        meta_dict = dict(refseq_accession = refseq_accession, 
                         collection_date = biosample.get('collection_date'), 
                         geo_loc_name = biosample.get('geo_loc_name'), 
                         isolate = biosample.get('isolate'), 
                         isolation_source = biosample.get('isolation_source'), 
                         host = biosample.get('host'), strain = biosample.get('strain'))
    else:
        print(f'failed to get json for {refseq_accession}, {status}')
    return meta_dict
        
def meta_dict_to_string(dict_, dict_keys):
    str_ = ''
    for i,j in enumerate(dict_keys):
        value_ = dict_[j]
        if value_ == None:
            value_ = 'NA'
        if i < len(dict_keys) -1:
            str_ += value_ + '\t'
        else:
            str_ += value_ + '\n'
    return str_

def download_multiple_refseq_meta(input_file, output_tsv, accession_idxs):
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
    #dict_list = []
    key_order = ['refseq_accession', 'collection_date', 'geo_loc_name', 'isolate', 'isolation_source', 'host', 'strain'] # look at the API page to get an idea of which fields to retrieve
    tsv_header = '\t'.join(key_order)
    f = open(output_tsv, 'a')
    f.write(tsv_header + '\n')
    for i,gcf in enumerate(gcfs):
        print(i,gcf)
        meta_dict = download_refseq_meta(gcf)
        if meta_dict != None:
            meta_str = meta_dict_to_string(meta_dict, key_order)
            #dict_list.append(meta_dict)
            f.write(meta_str)
    
    #df = pd.DataFrame.from_records(dict_list)
    #table = pd.read_table(meta_json)
    #df.to_csv(output_tsv, sep='\t', index=False)
    f.close()
def main(args):
    download_multiple_refseq_meta(args.input_file, args.output_tsv, args.accession_idxs)

if __name__ == '__main__':
    main(args)
