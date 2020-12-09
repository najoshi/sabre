> This is a fork of the [original repo](https://github.com/najoshi/sabre). I might be taking this tool into a different direction to what was originally intended

# sabre

> A cellular barcode demultiplexing tool of FASTQ files

## Content

- [Install](#install)
- [Quick start](#quick-start)
- [Usage](#usage)

## Install

```BASH
git clone https://github.com/serine/sabre
cd src
make
```

## Quick start

```BASH
sabre -f MultiplexRNASeq_S1_R1_001.fastq.gz \
      -r MultiplexRNASeq_S1_R2_001.fastq.gz \
      -b barcodes.txt \
      -c \
      -u \
      -m 2 \
      -l 10 \
      -a 1 \
      -s sabre.txt \
      -t 12
```

## Usage

> This tool is under development and this is very much an alpha version
> In it's current form the tool is highly customised a particular multiplexing protocol

### Cellular barcodes

In order to demultiplex the use needs to provide `barcodes.txt` file, which is three column tab delimited file

```
sample_name group barcode
```

currently group is semi-redundant column, it there for a feature that in the development. for most use cases group can equals to barcode

e.g

```
cntr_rep1    TAAGGCGA        TAAGGCGA
cntr_rep2    CGTACTAG        CGTACTAG
treat_rep1   AGGCAGAA        AGGCAGAA
treat_rep2   TCCTGAGC        TCCTGAGC
```
