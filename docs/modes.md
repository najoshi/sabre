# Sabre

## Different running modes

DOCS: In each case BARCODE and/or UMI are trimed off and
put into FASTQ header:

Not sure if I should have:

       BARCODE always has a precedent i.e BARCODE:UMI
       OR
       It follows the same structure as per experiment i.e
       if BARCODE+UMI then BARCODE:UMI
       else if UMI+BARCODE then UMI:BARCODE

All modes that begin with 3 will return single - R1 file, merging
R1 read into R2 header and renaming R2 into R1

10 = single-end where R1 has the following structure:

         R1 -->
         BARCODE+READ

20 = paired-end where R1 and R2 have the following structure:

         R1 -->                 <--R2
         BARCODE+READ----READ+BARCODE

this mode returns single file (R1) with barcode appended and into R1 header

30 = paired-end where R1 and R2 have the following structure:

         R1 -->     <-R2
         BARCODE----READ

40 = paired-end where

11 = single-end where R1 has the following structure:

         R1 -->
         BARCODE+UMI+READ

21 = paired-end where R1 and R2 have the following structure:

         R1 -->                         <--R2
         BARCODE+UMI+READ----READ+UMI+BARCODE

this mode returns single file (R1) with barcode appended and into R1 header

31 = paired-end where R1 and R2 have the following structure:

         R1 -->         <-R2
         BARCODE+UMI----READ

NOTE this gives me room for yet another mode e.g 12, 22, 32

12 = sinle-end where R1 has the following structure:

         R1 -->
         UMI+READ

22 = paired-end where R1 and R2 have the following structure:

         R1 -->         <--R2
         UMI+READ----READ+UMI

this mode returns single file (R1) with barcode appended and into R1 header

32 = paired-end where R1 and R2 have the following structure:

         R1 --> <-R2
         UMI----READ
