#!/bin/bash

# Organize variables.
column=$1 # column for which to calculate z-score
path_table_original=$2 # original table
path_table_novel=$3 # novel table
report=$4 # whether to print reports

# Set column delimiter.
# Calculate cumulative statistics.
count=$(cat $path_table_original | awk 'BEGIN { FS=" " } NR > 1 { n++ } END { print n }')
sum=$(cat $path_table_original | awk -v column=$column 'BEGIN { FS=" " } NR > 1 { sum += $column } END { print sum }')
sum_squares=$(cat $path_table_original | awk -v column=$column 'BEGIN { FS=" " } NR > 1 { sum_squares += ($column)^2 } END { print sum_squares }')
mean=$(cat $path_table_original | awk -v column=$column 'BEGIN { FS=" " } NR > 1 { sum += $column; n++ } END { if (n > 0) print (sum / n); }')
variations=$(cat $path_table_original | awk -v column="$column" -v mean="$mean" \
'BEGIN { FS=" " } NR > 1 { value += (($column - mean)^2) } END { print value }')
# The following formulas for standard deviation are equivalent.
#standard_deviation=$(echo "scale=3; sqrt(($sum_squares - (($sum^2)/$count))/($count - 1))" | bc -l)
#standard_deviation="echo 'sqrt(($sum_squares - (($sum^2)/$count))/($count - 1))' | bc -l"
#standard_deviation=$(bc -l <<< "sqrt(($variations)/($count - 1))")
standard_deviation=$(awk -v variations=$variations -v count=$count 'BEGIN {print sqrt(variations / (count - 1)) }')

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Note that this program currently is hard-coded for 6 column tables with"
  echo "the column to z-score in the 5th column. Not ideal."
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Initial statistics for column ${column} before z-score standardization..."
  echo "----------"
  echo "Count: ${count}"
  echo "Sum: ${sum}"
  echo "Sum of squares: ${sum_squares}"
  echo "Mean: ${mean}"
  echo "Variations: ${variations}"
  echo "Standard deviation: ${standard_deviation}"
fi

# Determine whether standard deviation is close to zero.
minimum=(-0.0001)
maximum=0.0001
zero_min=$(awk -v value=$standard_deviation -v minimum=$minimum \
'BEGIN {print (value>minimum) }')
zero_max=$(awk -v value=$standard_deviation -v maximum=$maximum \
'BEGIN {print (value<maximum) }')
if [[ $zero_min -eq 1 ]] && [[ $zero_max -eq 1 ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Standard deviation is close to zero."
  echo "Mean center values only."
  cat $path_table_original | awk 'BEGIN { FS=OFS=" " } NR == 1' > $path_table_novel
  cat $path_table_original | \
  awk -v column="$column" -v mean="$mean" -v \
  'BEGIN { FS=" "; OFS=" " } NR > 1 {print $1, $2, $3, $4, ($column - mean), $6}' >> $path_table_novel
else
  echo "----------"
  echo "----------"
  echo "----------"
  echo "Standard deviation is not close to zero."
  echo "Z-score standardize values."
  #head $path_table_original
  #rm $path_table_novel
  cat $path_table_original | awk 'BEGIN { FS=OFS=" "} NR == 1' > $path_table_novel
  cat $path_table_original | \
  awk -v column="$column" -v mean="$mean" -v stdev="$standard_deviation" \
  'BEGIN { FS=" "; OFS=" " } NR > 1 {print $1, $2, $3, $4, (($column - mean) / stdev), $6}' >> $path_table_novel
fi

# Set column delimiter.
# Calculate cumulative statistics.
count=$(cat $path_table_novel | awk 'BEGIN { FS=" " } NR > 1 { n++ } END { print n }')
sum=$(cat $path_table_novel | awk -v column=$column 'BEGIN { FS=" " } NR > 1 { sum += $column } END { print sum }')
sum_squares=$(cat $path_table_novel | awk -v column=$column 'BEGIN { FS=" " } NR > 1 { sum_squares += ($column)^2 } END { print sum_squares }')
mean=$(cat $path_table_novel | awk -v column=$column 'BEGIN { FS=" " } NR > 1 { sum += $column; n++ } END { if (n > 0) print (sum / n); }')
variations=$(cat $path_table_novel | awk -v column="$column" -v mean="$mean" \
'BEGIN { FS=" " } NR > 1 { value += (($column - mean)^2) } END { print value }')
# The following formulas for standard deviation are equivalent.
#standard_deviation=$(echo "scale=3; sqrt(($sum_squares - (($sum^2)/$count))/($count - 1))" | bc -l)
#standard_deviation="echo 'sqrt(($sum_squares - (($sum^2)/$count))/($count - 1))' | bc -l"
#standard_deviation=$(bc -l <<< "sqrt(($variations)/($count - 1))")
standard_deviation=$(awk -v variations=$variations -v count=$count 'BEGIN {print sqrt(variations / (count - 1)) }')

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Final statistics for column ${column} after mean center or z-score standardization..."
  echo "----------"
  echo "Count: ${count}"
  echo "Sum: ${sum}"
  echo "Sum of squares: ${sum_squares}"
  echo "Mean: ${mean}"
  echo "Variations: ${variations}"
  echo "Standard deviation: ${standard_deviation}"
  echo "----------"
  echo "----------"
  echo "----------"
fi
