#!/bin/bash

# need to remember to also remove variant count outliers

zcat $arrayvcf | grep -m 1 ^#CHROM | tr '\t' '\n' | tail -n +10 | sort > array.samples.sorted
cat $exchipvcf | grep -m 1 ^#CHROM | tr '\t' '\n' | tail -n +10 | sort > exchip.samples.sorted

cat $wgsmeta | tail -n +2 | awk '{print $4"-"$6}' | sort > wgs.samples.sorted
cat $wesmeta | tail -n +2 | awk '{print $4"-"$6}' | sort > wes.samples.sorted
cat $wesflag | cut -f1 | tail -n +3 | sort > wes.flagged.samples.sorted

comm -23 wes.samples.sorted wes.flagged.samples.sorted > keep1
comm -12 keep1 wgs.samples.sorted > keep2
wc -l keep2 # 139
comm -12 array.samples.sorted keep2 | wc -l # 0 
comm -12 exchip.samples.sorted array.samples.sorted | wc -l # 0
comm -12 exchip.samples.sorted keep2 | wc -l # 0
wc -l array.samples.sorted # 455
wc -l exchip.samples.sorted # 190 
