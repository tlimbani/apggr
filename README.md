# **Analysis of Population Genetics for Genomic Region (APGGR)**

A package with a set of modules for the analysis of population genetics of genomic regions targeting bacterial species

## Installation

pip install apggr

## Usage

Usage: apggr [--version\|--version-only] [--help] <command> <argument>

### **Commands**

cds

cluster

variant

### **Arguments**

--input-file

-d

--out

--out-file

--csv

## Examples

### BASH

apggr cds -f '/path/to/annotation/file' –csv 'path/to/csv'

### Python

from apggr.cds import get_cds_features
