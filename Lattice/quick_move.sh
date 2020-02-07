#! /bin/sh 

for file in *-0field.txt
  do
    echo "$file"
    filename=$file
    trimfile=$(echo "$filename" | cut -f 1 -d '.')
    newfile="${trimfile}-NoPot.txt"
    mv $filename $newfile
    echo "$newfile"
  done
