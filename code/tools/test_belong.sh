#####################################
# test_belong
# 
# This script checks whether a source file in a subdirectory has a matching file in model/src/
#
#------------------------------------


rac=$1;

echo "------------------------------"
echo " search for $rac"
echo " " 

echo " TEST 1: Does it appear in the model?"
  echo "           "

lsub=$(cat model/src/"$rac".f90 | grep -iv "END" | grep -i "subroutine" | sed "s/subroutine//g;s/SUBROUTINE//g" | cut -d"(" -f1)

for sub in $lsub
 do
  echo "           "
  echo "XXXXXXXX   " $sub
  grep -i "$sub" model/src/*
done

  echo "           "
echo " TEST 2:Does it appear in the list of objects?"
  echo "           "
cat optfile | grep "$rac"
