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

## Basic usage
presig is primarily designed to be used at the command line, though you can also import its individual functions
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

## Questions and bug reports
Please post an issue on the [GitHub](https://github.com/edawson/presig) if you have a question or find a bug.