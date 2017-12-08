############################
# monitor_expe.sh
#
# This script compares files that are found in an experiment subdirectory to the corresponding default model versions. 
# --------------------------
#
# V1.0 2013/03/05. JBG.

###########################


nb_arg=$#
exped=$1

# Number of arguments
if [ "$nb_arg" != "1" ]; then
  echo "Error: usage: sh tools/monitor_expe.sh subdirectory"
  exit
fi

# For all files that are found in the experiment directory,
lfic=$(ls $exped/??c/*)
for fic in $lfic
 do
  f1=$(echo $fic | sed "s/$exped/model/g")
  f2="$fic"
  echo " ";echo "####################################################################";echo "# $f2";echo "#-------------------";echo " ";
  sdiff $f1 $f2 -w 280 | colordiff
  echo " ";echo " ";echo " ";
done


