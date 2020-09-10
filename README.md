# ICTV Scripts
Scripts for parsing ICTV forms

# phage-sum

Use this script to extract GenBank (INSDC) summary information from the ICTV Viral Metadata Resource [VMR](https://talk.ictvonline.org/taxonomy/vmr/). It uses the accession number information to get metadata and genome information from the published sequence data, and returns a summary table.

## How to install

Clone the ictv-scripts repository into your folder of interest.
Install the necessary dependencies listed in the requirements.txt file.

```
pip install -r requirements.txt
```

## How to run
Here's the help:
```
usage: phagesum.py [-h] [-v] [--version] [--email EMAIL] [--csv]
                   [-o OUTPUT_FILE] [--accession_column ACCESSION_COLUMN]
                   summary_file

phagesum Summarise table of phage records and extract metrics from genome
annotation

positional arguments:
  summary_file          Intial summary file

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         verbose output
  --version             show program's version number and exit
  --email EMAIL         Email to query Entrez
  --csv                 Use this flag if your input is csv
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output file path
  --accession_column ACCESSION_COLUMN
                        Which column has the accession code

Licence: GPLv3 by Nabil-Fareed Alikhan <nabil@happykhan.com>
```
Your input file will be the a .txt or .csv (tabular) version of the VMR or a subset number of rows, including the header row. I don't recommend to use the whole VMR because that will take forever to run.
Additionally, you might want to remove the spaces in the accession column header (just to play it safe).

The columns of the output file are unorganised, sorry 'bout that.  
