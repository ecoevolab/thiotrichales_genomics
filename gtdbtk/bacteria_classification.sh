echo -e "genome_id\tclassification" > gtdbtk_classification.csv

for d in */; do
  file=$(ls "$d"/*.bac120.summary.tsv 2>/dev/null)
  if [ -f "$file" ]; then
    awk -F'\t' 'NR==2 {print $1 "\t" $2}' "$file"
  fi
done >> gtdbtk_classification.csv
