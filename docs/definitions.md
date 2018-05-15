define blocks along the read

BARCODE 
UMI
READ

then set values to different block

BARCODE = 8
UMI = 10

## Andrew's suggestion

```
--input sample_A_R1.fastq.gz:i8{index1},r151{read1},i8{index2}
```

```
--fq1 sample_A_R1.fastq.gz:i8{index1},r151{READ1},i8{index2}

--fq2 sample_A_R2.fastq.gz:i8{index1},r151{READ1},i8{index2}
```

We need to check that BARCODE == index1 in both fq1 and fq2 but also check that index1_fq1 == index1_fq2

```
--merge 12 merge R1 into R2
--merge 21 merge R2 into R1
```

either way resulting read is R1

```
--fq1 sample_A_R1.fastq.gz:8index1,*index2

--fq2 sample_A_R2.fastq.gz:i8{index1},r151{read2},i8{index2}
```
