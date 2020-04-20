presig: Transform a MAF file into the inputs of signet / sparsesigs
---------------------------------------
Eric T Dawson  
Nov 2019

![](https://github.com/edawson/presig/workflows/Presig%20Tests/badge.svg)


## Introduction
[presig](https://github.com/edawson/presig) generates count matrices of single-base (SBS) and
indel (ID) features for input to SigProfiler, SignatureAnalyzer, and signet/sparsesigs. It takes a
MAF file as input, as well as a corresponding references FASTA file.

## Requirements
presig requires pyfaidx for parsing FASTA files and pycotap for testing.

## Installation
```
git clone --recursive https://github.com/edawson/presig
cd presig/
pip install -r requirements.txt -e .
```

## Basic usage
presig is primarily designed to be used at the command line,
though you can also import its individual functions
into a python environment if desired.

### Convert a MAF file to feature counts
```
python presig/presig.py -m <MAF> -f <FASTA>
```

This will generate two files (`<MAF>.SBS96.tsv` and `<MAF>.ID83.tsv`) in the current directory.

### Convert a MAF file to the SigProfiler simple text format
```
python presig/presig.py -m <MAF> -f <FASTA> -s -u > <outputfile>.txt
```

### Run SigProfiler (including installing SigProfilerHelper)
To run SBS96 signatures using the TSV counts matrix generated using presig:
```
git clone --recursive https://github.com/edawson/sigprofilerhelper
python sigprofilerhelper/run_sigprofiler.py -c 16 -i 1000 -s 1 -e 7 -t <MAF>.SBS96.tsv
```

## Questions and bug reports
Please post an issue on the [GitHub](https://github.com/edawson/presig) if you have a question or find a bug.
