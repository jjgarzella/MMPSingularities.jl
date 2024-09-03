#!/bin/bash

FILE="heightsbargraph.csv"
awk -F, '
NR==1 { header = $0; next }
{
  for (i=1; i<=NF; i++) sum[i] += $i
}
END {
  print header
  for (i=1; i<=length(sum); i++) {
    printf sum[i]
    if (i < length(sum)) printf ","
  }
  printf "\n"
}' $FILE > tmpfile && mv tmpfile $FILE
