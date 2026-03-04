#!/usr/bin/env bash
set -euo pipefail

module purge
module load miller/6.16.0

BASEDIR="$(dirname "$0")"

# Must provide at least one library ID
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 LIBRARY_ID [LIBRARY_ID ...]" >&2
    exit 1
fi

# Get latest catalogs
FASTQS_CATALOG=(/datasets/caeg_production/_STATS/20*.fastq.tsv)
RESULTS_CATALOG=(/datasets/caeg_production/_STATS/20*.prod.tsv)

FASTQS_CATALOG="$(realpath "${FASTQS_CATALOG[-1]}")"
RESULTS_CATALOG="$(realpath "${RESULTS_CATALOG[-1]}")"

# Temporary file containing library list
SEARCH_LIST=$(mktemp)

{
    echo "library"
    for lib in "$@"; do
        echo "$lib"
    done
} > "$SEARCH_LIST"

# Join and filter
mlr --tsv --from "$FASTQS_CATALOG" \
    join -f "$RESULTS_CATALOG" -j library,date,flowcell --lp results_ --rp fastq_ \
    then join -j library -f "$SEARCH_LIST"

rm "$SEARCH_LIST"