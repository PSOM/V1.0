
mv -f _check.tmp poub

lfic=$(ls log_Z*)

for fic in $lfic
 do
  o1=$(ls -al $fic)
  o2=$(cat $fic | grep "kin" | tail -1)
  echo $o2 "  -  " $o1 >> _check.tmp
done

cat _check.tmp | sort
