Calculation of pop gen stats was performed in Stacks's populations (v1.44) module using the filtered freebayes vcf.

FST

Calculation of window -based FST and assessment of significance followed the approach of Hohenlohe (2010; https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000862) and the Stacks manual. A challenging aspect of RAD data is that the number of SNPs in a given genomic window can vary quite a lot and therefore assessment of significance should probably take that into account. The Hohenlohe approach performs bootstrap resampling to test the significance of a given window smoothed FST against an empirical distribution that reflects the characteristics of the window and the genomewide FST. Given that bootstrap resampling is computationally expensive, it is recommended to perform a low number of replicates initially and then increase the number of replicates for windows at the tails of the distribution as necessary. We initially performed 100,000 replicates for all windows and then created a "whitelist" of windows (p < 0.001) to perform the bootstrap resampling with 10 million replicates. 

Run with 100k replicates:

`
populations -t 1 -b 1 -V SNP.msp.recode.vcf --out_path ./100k/ -M MSP.popmap -f p_value -r 0.75 --bootstrap_fst --fstats --bootstrap_reps 100000 -k --window_size 25000
`

Create a whitelist for sites with pvalue < 0.001 in the initial output:

`
cat SNP.msp.recode.fst_INS-OUT.tsv | awk '$20 < 0.001' | cut -f2 > white.list
`

Run with 10M replicates:

`
populations -t 1 -b 1 -V SNP.msp.recode.vcf --out_path ./100k/10m/ -M MSP.popmap -f p_value -r 0.75 --bootstrap_fst --fstats --bootstrap_reps 10000000 -k --window_size 25000 --bootstrap_wl ./100k/white.list
`

Nucleotide diversity is also calculated in Stacks populations. See Rcode for Z-Score-based analysis of nucleotide diversity rations looking for population specific loss of diversity. NOTE: These are variant site nucleotide diversity estimates (i.e. does not account for adjacent invariant sites). 
