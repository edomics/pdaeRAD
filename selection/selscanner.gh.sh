#!/bin/bash

#Create temporary vcfs for each scaffold and loop over them to calculate ehh and xpehh
while read line;
do
grep $line ../SNP.dom.phased.INS.recode.vcf > temp.ins.vcf
grep $line ../SNP.dom.phased.OUT.recode.vcf > temp.out.vcf
grep $line ../SNP.msp.recode.phased.bi.map > temp.map

selscan --ihs --vcf temp.ins.vcf --map temp.map --out $line"ins" --trunc-ok
selscan --ihs --vcf temp.out.vcf --map temp.map --out $line"out" --trunc-ok
selscan --xpehh --vcf temp.ins.vcf --vcf-ref temp.out.vcf --map temp.map --out $line"ivso" --trunc-ok
done < x20snps.300k.scafs


