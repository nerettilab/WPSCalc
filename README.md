# WPSCalc
Calculates Window-Protection Score from aligned sequencing data over a given Chromosome &amp; interval

## Instructions:
Specify bedfile, reference genome, chromosome, start coordinate, end coordinate, window size, and paired-end (True/False) as input arguments.

For example:
> python WPSCalc.py test.bed hg38 6 1000000 1250000 50 --pairedend True

...will run the script on a paired-end bedfile called 'test.bed' using hg38, returning scores for chr6:1000000-1250000 with a window size of 50

The script will output a list of basepairs and corresponding scores over the specified range.
