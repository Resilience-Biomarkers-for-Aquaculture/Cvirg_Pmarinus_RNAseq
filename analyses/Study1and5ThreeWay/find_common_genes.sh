#!/usr/bin/env bash
# Usage: ./common_gene_ids.sh file1.tsv file2.tsv ...

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 file1.tsv file2.tsv ..." >&2
    exit 1
fi

# Start with gene_ids from the first file
cut -f1 "$1" | sort -u > tmp_common_ids

# Intersect with gene_ids from the remaining files
shift
for file in "$@"; do
    cut -f1 "$file" | sort -u > tmp_ids
    comm -12 tmp_common_ids tmp_ids > tmp_common_new
    mv tmp_common_new tmp_common_ids
done

# Output the result
cat tmp_common_ids

# Cleanup
rm tmp_common_ids tmp_ids 2>/dev/null
