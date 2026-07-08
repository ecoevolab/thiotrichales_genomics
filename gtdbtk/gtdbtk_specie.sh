awk '
/^[0-9]/ {
  if (NR>1) print record
  record=$0
  next
}
{
  record=record " " $0
}
END {
  print record
}
' gtdbtk_classification.csv \
| awk '
{
  id=$1
  if (match($0, /s__[^;]+/, a)) {
    sub(/^s__/, "", a[0])
    print id "\t" a[0]
  }
}
' \
| sort -k1,1n > genome_species_sorted.csv
