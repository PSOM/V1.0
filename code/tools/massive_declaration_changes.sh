
dir1=$(cat ./optfile | grep -v "#" | grep "dirc=" | cut -d= -f2)

#dir=$dir1"/code/model/src/"
dir=$dir1"/code/model/inc/"

lfic=$(ls $dir/*.f90)

for fic in $lfic
 do
  echo $fic
  cat $fic | sed "s/REAL(/@1/g;s/REAL,/@2/g;s/really/@3/g" | sed "s/REAL/real/g" | sed "s/ real\*8[ ]*::/ @41 /g;s/ real\*4[ ]*::/ @41 /g;s/ real[ ]*::/ @41 /g;s/ real\*8/ @41 /g;s/ real\*4/ @41 /g;s/ real / @41 /g" | sed "s/double precision /@41 /g" | sed "s/\@41/REAL(kind=rc_kind) ::/g" | sed "s/@1/REAL(/g;s/@2/REAL,/g;s/@3/really/g" > _fic.tmp
  sleep 1
  cat _fic.tmp | sed "s/ :: ,/,/g" > _fic2.tmp
  sleep 1
  cp _fic2.tmp $fic

#  cat $fic | sed "s/double precision/REAL*8 :: /g"
#cat $fic
done
