#!/bin/bash
# Compare CPU miRanda vs CUDAMiranda on the same inputs.
M=/home/user/miniconda3/envs/miranda_env/bin/miranda
G=/mnt/e/work_KI/cudaMiranda/compat/CUDAMiranda
D=/mnt/e/work_KI/cudaMiranda/data
O=/tmp/cmp; mkdir -p $O
for t in "query_1 ref_1000_256" "query_32 ref_1000_256" "query_256 ref_1000_256"; do
  set -- $t; n=$1_$2
  s=$(date +%s.%N); $M $D/$1.fa $D/$2.fa -out $O/$n.cpu >/dev/null 2>&1; c=$(echo "$(date +%s.%N)-$s"|bc)
  s=$(date +%s.%N); $G $D/$1.fa $D/$2.fa -out $O/$n.gpu >/dev/null 2>&1; g=$(echo "$(date +%s.%N)-$s"|bc)
  # per-hit lines (single '>' followed by a name), and per-pair summaries ('>>')
  for k in cpu gpu; do
    grep -P '^>(?!>)' $O/$n.$k | sort > $O/$n.$k.hits
    grep '^>>' $O/$n.$k | sort > $O/$n.$k.sum
  done
  ch=$(wc -l < $O/$n.cpu.hits); gh=$(wc -l < $O/$n.gpu.hits)
  common=$(comm -12 $O/$n.cpu.hits $O/$n.gpu.hits | wc -l)
  sc=$(comm -12 $O/$n.cpu.sum $O/$n.gpu.sum | wc -l)
  printf "%-24s cpu %7.1fs gpu %6.1fs | hits cpu=%d gpu=%d identical=%d | pair summaries cpu=%d gpu=%d identical=%d\n" \
    $n $c $g $ch $gh $common $(wc -l < $O/$n.cpu.sum) $(wc -l < $O/$n.gpu.sum) $sc
done
