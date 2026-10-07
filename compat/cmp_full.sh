#!/bin/bash
M=/home/user/miniconda3/envs/miranda_env/bin/miranda
G=/mnt/e/work_KI/cudaMiranda/compat/CUDAMiranda
D=/mnt/e/work_KI/cudaMiranda/data
O=/tmp/cmp; mkdir -p $O
for q in query_1 query_32; do
  n=${q}_utr
  s=$(date +%s); $G $D/$q.fa $D/utr3prime_targetscan.fa -out $O/$n.gpu -quiet >/dev/null 2>$O/$n.gpu.err; ge=$?; g=$(( $(date +%s)-s ))
  s=$(date +%s); $M $D/$q.fa $D/utr3prime_targetscan.fa -out $O/$n.cpu -quiet >/dev/null 2>&1; c=$(( $(date +%s)-s ))
  for k in cpu gpu; do grep -P '^>(?!>)' $O/$n.$k | sort > $O/$n.$k.hits; grep '^>>' $O/$n.$k | cut -f1-6,8- | sort > $O/$n.$k.sum; done
  echo "$n: gpu exit=$ge ${g}s, cpu ${c}s | hits cpu=$(wc -l <$O/$n.cpu.hits) gpu=$(wc -l <$O/$n.gpu.hits) identical=$(comm -12 $O/$n.cpu.hits $O/$n.gpu.hits|wc -l) | summaries identical=$(comm -12 $O/$n.cpu.sum $O/$n.gpu.sum|wc -l)/$(wc -l <$O/$n.cpu.sum)"
  cat $O/$n.gpu.err
done
