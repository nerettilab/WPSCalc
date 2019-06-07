# WPSCalc
Calculates Window-Protection Score (WPS) from aligned sequencing data over a given chromosome &amp; interval

WPS is defined for each basepair within a region, for a given window of size K. The score for a given base is defined as the number of molecules spanning the full window (basepair ± K/2) minus the number of molecules that align with an endpoint inside the window.

## Instructions:
Specify bedfile, reference genome, chromosome, start coordinate, end coordinate, window size, and paired-end (True/False) as input arguments.

For example:
> python WPSCalc.py test.bed hg38 6 1000000 1250000 50 --pairedend True

...will run the script on a paired-end bedfile called 'test.bed' using hg38, returning scores for chr6:1000000-1250000 with a window size of 50

The script will output a list of basepairs and corresponding scores over the specified range.
