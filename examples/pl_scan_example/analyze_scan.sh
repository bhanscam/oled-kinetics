#!/bin/bash

# bash analyze_scan.sh N_nodes out_file_name

echo "" > $2
i=0

for ((node=1; node<=$1; node++))
do
for core in {1..128}
do

if test -f run-$node-$core.log
then

if ! grep -q "too big to handle" run-$node-$core.log
then
j=$(tail -1 run-$node-$core.log)
echo $j >> $2
fi

i=$[$i+1]
fi

done
done

sort -o  $2 $2
echo "total run files found: " $i
