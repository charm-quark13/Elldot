#! /bin/sh 

for file in discrete_*.o*
  do
    echo "Calculating $file"
    ./$file
    echo "$file completed"
  done
