#####################
# tommy.sh
#
# Little tool that allows to create sensitivity experiments using the "user" namelist.
#
# JBG. 03/06/2013. 
#
#-------------------


# Initialization

root= #Experiment name
proc_nb=4
dirc=$(cat optfile | sed "s/^dirc=/@=/g" | grep "@=" | cut -d= -f2)

cat > f_tommy.tmp << EOF

#standard
1.  1. 1. 1. 2000. 1. 1. 1., 
         
# position of wiggle changed
1.  1. 1. 1. 2005. 1. 1. 1., 
1.  1. 1. 1. 2015. 1. 1. 1., 
1.  1. 1. 1. 2025. 1. 1. 1., 
1.  1. 1. 1. 1995. 1. 1. 1., 
1.  1. 1. 1. 1985. 1. 1. 1., 
                
# amplitude of wiggle changed
1.  1. 1. 1. 2000. 0.5 1. 1., 
1.  1. 1. 1. 2000. 0.2 1. 1., 
1.  1. 1. 1. 2000. 0. 1. 1., 

EOF


# Initialization is over.
#------------------------


lconfs=$(cat f_tommy.tmp | grep -v "#")
lconfs=$(echo $lconfs | sed "s/,[ ]*/X/g" | sed "s/ /_/g" | sed "s/X/ /g")
i=0;n=1;

datet=$(date | sed "s?  ? 0?g" | sed "s? ?_?g" | awk -F_ '{print $6"_"$2"_"$3"_"$4}' | sed "s?:?_?g")
lexec="./lexec_"$datet"__"$n
lexec_pr=$(ls lexec* | wc -w)

# Checks that no previous script files are still in ./code
if [ $lexec_pr = 0 ];
 then
  echo "OK"
  else
   echo $lexec_pr 
  exit
fi

# End of pre-processing
#----------------------


for conf in $lconfs   # For every sensitivity experiment,
 do

  let i=$i+1

  u1=$(echo $conf | cut -d_ -f1)
  u2=$(echo $conf | cut -d_ -f2)
  u3=$(echo $conf | cut -d_ -f3)
  u4=$(echo $conf | cut -d_ -f4)
  u5=$(echo $conf | cut -d_ -f5)
  u6=$(echo $conf | cut -d_ -f6)
  u7=$(echo $conf | cut -d_ -f7)
  u8=$(echo $conf | cut -d_ -f8)

  u1p=$(echo $u1  | sed "s/\-/m/g" | sed "s/\./p/g")
  u2p=$(echo $u2  | sed "s/\-/m/g" | sed "s/\./p/g")
  u3p=$(echo $u3  | sed "s/\-/m/g" | sed "s/\./p/g")
  u4p=$(echo $u4  | sed "s/\-/m/g" | sed "s/\./p/g")
  u5p=$(echo $u5  | sed "s/\-/m/g" | sed "s/\./p/g")
  u6p=$(echo $u6  | sed "s/\-/m/g" | sed "s/\./p/g")
  u7p=$(echo $u7  | sed "s/\-/m/g" | sed "s/\./p/g")
  u8p=$(echo $u8  | sed "s/\-/m/g" | sed "s/\./p/g")

  key=$root"_"$u1p"_"$u2p"_"$u3p"_"$u4p"_"$u5p"_"$u6p"_"$u7p"_"$u8p


  # Creation of a new namelist
  nm_name="namelist_"$key
  cp namelist_"$root" $nm_name 

  # Creation of a new output directory
  dirout="$dirc/output_"$key

  echo "-- "
  echo "namelist: " $nm_name
  echo "dirout  : " $dirout

  mkdir $dirout > /dev/null 2> /dev/null

  # Checks whether the output directory is emply
  files_ls=$(ls $dirout)

  files_nb=$(echo $files_ls | wc -w)
  if [ $files_nb = 0 ]; then
    res=$(echo "")
   else
    res=$(echo ": DANGER")
  fi

  # Modification of the standard namelist into the experiment namelist
  echo "Creation of dirout. "$files_nb "file(s)" $res
  cat "namelist_"$root | sed "s%$root%$key%g" > $nm_name
  cat $nm_name | sed "s/user1=X/user1=$u1/g" > $nm_name".tmp"
  cat $nm_name".tmp" | sed "s/user2=X/user2=$u2/g" > $nm_name
  cat $nm_name | sed "s/user3=X/user3=$u3/g" > $nm_name".tmp"
  cat $nm_name".tmp" | sed "s/user4=X/user4=$u4/g" > $nm_name
  cat $nm_name | sed "s/user5=X/user5=$u5/g" > $nm_name".tmp"
  cat $nm_name".tmp" | sed "s/user6=X/user6=$u6/g" > $nm_name
  cat $nm_name | sed "s/user7=X/user7=$u7/g" > $nm_name".tmp"
  cat $nm_name".tmp" | sed "s/user8=X/user8=$u8/g" > $nm_name

  echo ""
  cat $nm_name | grep "user" | sed "s/^/| /g"
  echo ""

  # Writing of the script files that will contain lines like
  # ./exe/nh_expe < namelist_expe_1 > log_expe_1 &

  exe="./exe/nh_$root"
  code_ln=$(echo $exe "< $nm_name > log_$key" )

  echo $i; echo $code_ln;

  # Management of the number of experiment per file (given the number of processor of the machine)
  let proc_ab=$proc_nb+1
  if [ $i = $proc_ab ]; then
    let i=1
    let n=$n+1
    datet=$(date | sed "s?  ? 0?g" | sed "s? ?_?g" | awk -F_ '{print $6"_"$2"_"$3"_"$4}' | sed "s?:?_?g")
    echo "DD " $datet
    lexec="./lexec_"$datet"__"$n
  fi
   
  if [ $i = $proc_nb ]; then
    echo $code_ln | sed "s/$/ /g" >> $lexec
   else
    echo $code_ln | sed "s/$/ \&/g" >> $lexec
  fi

done # done for every experiment

#---------------------------------

#---------------------------------
# Final

lexec_ls=$(ls lexec*)
nexec=$(ls lexec* | wc -l)
chmod 777 lexec*
echo $nexec

lfinal=$(echo $lexec_ls | sed "s/ / ; /g")

echo " "
echo " "
echo "execution of $lfinal "

