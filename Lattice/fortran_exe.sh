#! /bin/sh 

for file in *.o*
  do
    echo "Calculating $file"
    ./$file
    echo "$file completed"
  done
