#!/bin/bash

#inputs = batchnum, n_mig, prefix
 
random=$(shuf -i 1-10000000 -n 1)

for i in $(eval echo "{0..$2}"); do  treemix -i ../$3.gz -m $i -o $3.noboot.m.$i.r.$1 -root TAB -seed $random > $3.r.$1.treemix_${i}_log
done


