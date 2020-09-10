#!/usr/bin/env python3
"""
phagesum Summarise table of phage records and extract metrics from genome annotation 

phagesum Summarise table of phage records and extract metrics from genome annotation 

### CHANGE LOG ### 
2020-09-10 Nabil-Fareed Alikhan <nabil@happykhan.com>
    * Initial build - split from dirty scripts
"""
import csv
from Bio import Entrez
from Bio import SeqIO
import pprint
import logging
import csv 
import datetime
from os import path, mkdir, listdir
import os 
import argparse
import meta
import sys
import time

epi = "Licence: " + meta.__licence__ +  " by " +meta.__author__ + " <" +meta.__author_email__ + ">"
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s')
log = logging.getLogger(__name__)

def main(args):
    Entrez.email = args.email
    new_dict = []

    with open(args.summary_file) as f:
        dia = csv.excel_tab
        if args.csv:
            dia = csv.excel
        file_list = csv.DictReader(f, dialect=dia)
        headers = []
        for line in file_list:
            new_line = line
            acc = line[args.accession_column]
            if acc == "":
                log.error('No accession code for ' + line.get('Virus name(s)', 'VIRUS NAME NOT FOUND') )
                continue
            with Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text") as genbank_handle:
                seq_record = SeqIO.read(genbank_handle, "genbank")
                source = [x for x in seq_record.features if x.type =='source'][0]
                for k,v  in dict(source.qualifiers).items():
                    if k == 'lab_host':
                        new_line['host'] = ','.join(v)
                    elif k == 'isolate':
                        new_line['strain'] = ','.join(v)
                    else:
                        new_line[k] = ','.join(v)
                new_line['topology'] = seq_record.annotations['topology']
                new_line['submission_date'] = seq_record.annotations['date']
                new_line['organism'] = seq_record.annotations['organism']
                new_line['trna_count'] = len([x for x in seq_record.features if x.type == 'tRNA'])
                new_line['cds_count'] = len([x for x in seq_record.features if x.type == 'CDS'])
                repeat_regions = [x for x in seq_record.features if x.type == 'repeat_region']
                for x in  repeat_regions:
                    if x.location.nofuzzy_start == 0:
                        new_line['terminal_repeat'] = 'Yes'

                new_line['repeat_region_count'] = len(repeat_regions)
                new_line['intron_count'] = len([x for x in seq_record.features if x.type == 'intron'])
                new_line['feature_types'] = ','.join(sorted(set([x.type for x in seq_record.features])))
                new_line['length'] = len(seq_record.seq)
                pmid_list = [x.pubmed_id for x in seq_record.annotations['references'] if x.pubmed_id != '']
                new_line['pmid'] = ','.join(pmid_list)
            headers += new_line.keys()
            headers = list(set(headers))
            new_dict.append(new_line)
            pprint.pprint(dict(new_line))
        file_out = csv.DictWriter(open(args.output_file, 'w'), fieldnames=headers)
        file_out.writeheader()
        file_out.writerows(new_dict)

def is_valid_file(parser, arg):
    if not path.exists(arg):
        parser.error("The File %s does not exist!" % arg)
    else:
        return arg


if __name__ == '__main__':
    start_time = time.time()
    log.setLevel(logging.INFO)
    desc = __doc__.split('\n\n')[1].strip()
    parser = argparse.ArgumentParser(description=desc,epilog=epi)
    parser.add_argument ('-v', '--verbose', action='store_true', default=False, help='verbose output')
    parser.add_argument('--version', action='version', version='%(prog)s ' + meta.__version__)
    parser.add_argument('summary_file', action='store', help='Intial summary file', type=lambda x: is_valid_file(parser, x))
    parser.add_argument('--email', action='store',  help='Email to query Entrez', default="nabil@happykhan.com")
    parser.add_argument('--csv', action='store_true',  help='Use this flag if your input is csv', default=False)
    parser.add_argument('-o', '--output_file', action='store',  help='Output file path', default="phagesum.csv")
    parser.add_argument('--accession_column', action='store',  help='Which column has the accession code', default="Virus GENBANK accession")
    args = parser.parse_args()
    if args.verbose: 
        log.setLevel(logging.DEBUG)
        log.debug( "Executing @ %s\n"  %time.asctime())    
    main(args)
    if args.verbose: 
        log.debug("Ended @ %s\n"  %time.asctime())
        log.debug('total time in minutes: %d\n' %((time.time() - start_time) / 60.0))
    sys.exit(0)