#!/bin/bash
# Run CUDAMiranda on a query FASTA of any size by splitting it into chunks of
# at most 480 miRNAs (the per-run limit) and concatenating the results.
# Usage: run_split.sh queries.fa targets.fa out.txt [extra CUDAMiranda options...]
set -e
[ $# -lt 3 ] && { echo "usage: $0 queries.fa targets.fa out.txt [options]"; exit 1; }
Q=$1; R=$2; OUT=$3; shift 3
BIN="$(cd "$(dirname "$0")" && pwd)/CUDAMiranda"
# Batch 16 references per kernel launch (~3-4x faster than 1). With 480 queries
# and 22.5 kb UTRs, -rbs 32 exceeds 8 GB of GPU memory; pass -rbs N to override.
[[ " $* " == *" -rbs "* ]] || set -- -rbs 16 "$@"
TMP=$(mktemp -d)
awk -v dir="$TMP" '/^>/{ if (n % 480 == 0) f = sprintf("%s/chunk_%04d.fa", dir, n / 480); n++ } { print > f }' "$Q"
: > "$OUT"
for c in "$TMP"/chunk_*.fa; do
  echo "Running $(basename "$c") ($(grep -c '^>' "$c") queries)" >&2
  "$BIN" "$c" "$R" -out "$c.out" -quiet "$@" > /dev/null
  cat "$c.out" >> "$OUT"
done
rm -rf "$TMP"
