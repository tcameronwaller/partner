#!/bin/bash

# Organize variables.
value=$1 # floating-point value for comparison

minimum=(-0.0001)
maximum=0.0001

# Compare value to minimal and maximal thresholds near zero.
#zero_min=$(awk -v value=$standard_deviation -v minimum=$minimum \
#'BEGIN {print value>minimum?1:0 }')
#zero_max=$(awk -v value=$standard_deviation -v maximum=$maximum \
#'BEGIN {print value<maximum?1:0 }')
zero_min=$(awk -v value=$value -v minimum=$minimum \
'BEGIN {print (value>minimum) }')
zero_max=$(awk -v value=$value -v maximum=$maximum \
'BEGIN {print (value<maximum) }')
echo "zero_min: " $zero_min
echo "zero_max: " $zero_max
if [[ $zero_min -eq 1 ]] && [[ $zero_max -eq 1 ]]; then
  echo "value ${value} is close to zero"
  echo "(between ${minimum} and ${maximum})"
else
  echo "value ${value} is not close to zero"
  echo "(not between ${minimum} and ${maximum})"
fi
