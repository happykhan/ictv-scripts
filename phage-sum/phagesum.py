import csv
from Bio import Entrez
from Bio import SeqIO
import pprint

Entrez.email = "nabil@happykhan.com"
new_dict = []

with open('/usr/users/QIB_fr005/alikhan/Downloads/Ino-Plectroviridae.csv') as f:
    file_list = csv.DictReader(f)
    headers = []
    for line in file_list:
        new_line = line
        acc = line['Accession']
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
    file_out = csv.DictWriter(open('/usr/users/QIB_fr005/alikhan/Downloads/Ino-Plectroviridae.done.csv', 'w'), fieldnames=headers)
    file_out.writeheader()
    file_out.writerows(new_dict)

