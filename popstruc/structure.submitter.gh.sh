#!/bin/bash

#inputs = kval startrep endrep

seq $2 1 $3 > index
n=$(wc -l index | awk '{print $1}') 
shuf -i 1-10000000 -n $n | paste index - > index.random

while read index random
do
echo "structure -i SNP.distLD.bi.struc.linux -o SNP.distLD.bi.$1.K.rep.$index -K $1 -D $random -m mainparams -e extraparams" | cat script.header - > structure.run.$1.K.rep.$index.sh
sbatch structure.run.$1.K.rep.$index.sh
done < index.random

