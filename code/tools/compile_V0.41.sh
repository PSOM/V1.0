############################
# compile_complex.sh
#
# This script enables the combined use of model and model_build 
#  to create a hybrid executable with both old and new files.
# --------------------------
#
# V1.0 2012/02/27. JBG.
#  
#




###########################



list_of_files2() {

echo "-------------------------------------------------------"
echo "ls src :"
echo "`ls -al src`"

echo "ls inc :"
echo "`ls -al inc`"

echo "-------------------------------------------------------"

}

list_of_files() {
echo ""
}

######################
# Arguments

nb_arg=$#
rep_alt=$1

if [ "$nb_arg" != "0" ] && [ "$nb_arg" != "1" ]; then
  echo "error: usage: sh compile*.sh directory"
  exit
fi

if [ "$nb_arg" = "0" ]; then
  echo "default setting."
  dir_alt="model_build"
  exe_ext=""
 else
  echo "custom setting."
  dir_alt=$rep_alt
  exe_ext="_"$dir_alt
fi


  echo "alternative directory       : " $dir_alt 



######################
# Reading of the optfile

optfile="./optfile"
dirc=$(cat $optfile | grep -v "#" | grep "dirc=" | cut -d= -f2) 

dirloc=`pwd`
EXE="nh"
exef="$EXE""$exe_ext"

echo "creation of the executable  : "$exef 

################################################

echo "###########################"
echo "# compile_complex.sh starts"
echo "# "
echo "# Use of model and model_build to create $exef"
echo "# "

dirh=$dirc"/code/model"

# Default directories model/src and inc
dirsrc=$(echo $dirh"/src")
dirinc=$(echo $dirh"/inc")

# Alternative directories model_build/src and inc
dirsrcalt=$(echo $dirh"/../$dir_alt/src")
dirincalt=$(echo $dirh"/../$dir_alt/inc")


# END OF THE HEADER
################################################
# START OF THE SCRIPT

##############
# Step1 : Detecting alternative files.

# List the alternative directories
lsrcalt=$(ls $dirsrcalt)
lincalt=$(ls $dirincalt)

list_of_files

echo "-------------------------------------------------------"
echo "These files are detected in the model_build directories:"
echo "/src : " $lsrcalt
echo "/inc : " $lincalt
echo "-------------------------------------------------------"

# Option: the presence of a file in one of the alternative directories, even of an old one is a new event.
# In order to match the "makefile" philosophy (based on dates), every file found in the alternative directories are touched.
for fic in $lsrcalt
 do
  touch $dirsrcalt/$fic
done

for fic in $lincalt
 do
  touch $dirincalt/$fic
done


# What follows is a security test:
# For every file that is found in the alternative src directory, we check its presence in the list of objects in ./optfile.

for srcalt in $lsrcalt
 do
  srcalto=$(echo "@"$srcalt"@" | sed "s/\.f90//g")
  srcalto2=$(echo $srcalto | sed "s/@//g")
  sloc=$(cat ./optfile | grep "="  | sed "s/=/=@/g;s/$/@/g" | sed "s/ /@/g" | sed "s/lobj/X@@X/g;s/lmod/X@@X/g" | sed "s/\.mod//g;s/\.o//g" | grep "X@@X" | cut -d= -f2 )
  sloc=$(echo $sloc | grep $srcalto)
  if [ "$sloc" = "" ]; then
    echo "Warning ! " $srcalto2 " does not appears in the list of objects nor in the list of modules in optfile. $srcalto2 will not be taken into account during compilation !!"
   else
    echo $srcalto2 "appears in the list of objects or modules in optfile. It will be taken into account during compilation."
  fi

done 





#############
# Step2 : Replacing temporarily the old files, if exist.

echo "-------------------------------------------------------"
echo "/src :"

for src in $lsrcalt
 do 
  # For every alternative file, the presence of the old file is checked.
  if [ -s $dirsrc/$src ]; then
     # If the old file exists, its name is temporarily changed
    echo "   $dirsrc/$src is moved into $dirsrc/$src""_oldccx"
    mv $dirsrc/$src $dirsrc/$src"_oldccx"
   else
    # If it does not exists, a file is created to remember that an alternative file is present.
    echo "Warning ! The file " $src " that has been found in " $dirsrcalt " has no equivalent in " $dirsrc
    echo "  Therefore, $src will simply be linked in $dirsrc . A file" $dirsrc/$src"_new_file" will be created to remember it.
    touch $dirsrc/$src"_new_file"
  fi
    # In every case, a link is made between model/src and model_build/src
    ln -s $dirsrcalt/$src $dirsrc/$src 
    echo "   /`ls -al $dirsrc/$src | cut -d/ -f2-`"
done

#Idem for include files
echo "/inc :"
for inc in $lincalt
 do  
  if [ -s $dirinc/$inc ]; then
    echo "   $dirinc/$inc is moved into $dirinc/$inc""_oldccx"
    mv $dirinc/$inc $dirinc/$inc"_oldccx"
   else
    echo "Warning ! The file " $inc " that has been found in " $dirincalt " has no counterpart in " $dirinc
    echo " A file" $dirinc/$src"_new_file" will be created to remember it.
    touch $dirinc/$inc"_new_file"
  fi  
    ln -s $dirincalt/$inc $dirinc/$inc 
    echo "   /`ls -al $dirinc/$inc | cut -d/ -f2-`"
done

list_of_files

##############
# Step3 : Compiling

echo "-------------------------------------------------------"
echo "Creation of the makefile "

$dirc/code/tools/genmakefilel

echo "Building the executable"

make -C $dirc/code/mkfile $EXE
#make -C $dirc/code/mkfile clean
cp $dirc/code/mkfile/$EXE $dirc/code/exe/$exef

echo "-------------------------------------------------------"

#############
# Step4 : Removing the temporary links, replacing the old files.

# 1. Standard case: both alternative and old files existed. 
# The _oldccx files are searched for.

lsrcmv=$(ls $dirsrc | grep "_oldccx")
lincmv=$(ls $dirinc | grep "_oldccx")

echo "-------------------------------------------------------"
echo " These files have been displaced :"
echo "/src : " $lsrcmv
echo "/inc : " $lincmv
echo "-------------------------------------------------------"

for srcmv in $lsrcmv
 do
  # For every files couple, the link is broken and the old file is named as before.
  src=$(echo $srcmv | sed "s/_oldccx//")
  rm -f $dirsrc/$src
  echo "$dirsrc/$src is back"
  mv $dirsrc/$srcmv $dirsrc/$src
done

# Idem for include files.
for incmv in $lincmv
 do
  inc=$(echo $incmv | sed "s/_oldccx//")
  rm -f $dirinc/$inc
  echo "$dirinc/$inc is back"
  mv $dirinc/$incmv $dirinc/$inc
done

# 2. Non-standard case: only the alternative file existed.
# The _new_file files are searched for. 

lsrcnewf=$(ls $dirsrc | grep "_new_file")
lincnewf=$(ls $dirinc | grep "_new_file")

lsrcnewf2=$(echo $lsrcnewf | sed "s/_new_file//g")
lincnewf2=$(echo $lincnewf | sed "s/_new_file//g")


echo "-------------------------------------------------------"
echo " Follows a list of files of ./model_build from which a copy has been made in ./model. These copies will now be removed." 
echo "/src : " $lsrcnewf2
echo "/inc : " $lincnewf2
echo "-------------------------------------------------------"

# The link is broken and the remainder is deleted. 
for fic in $lsrcnewf2
 do
  echo "The files " $dirsrc/$fic " and " $dirsrc/$fic"_new_file have been removed."
  rm -f $dirsrc/$fic
  rm -f "$dirsrc/$fic"_new_file
done

for fic in $lincnewf2
 do
  echo "The files " $dirinc/$fic " and " $dirinc/$fic"_new_file have been removed."
  rm -f $dirinc/$fic
  rm -f "$dirinc/$fic"_new_file
done


###############
#Step5: Ending 


list_of_files

echo "# "
echo "# `ls -al $dirc/code/exe/$exef`"
#echo "# $dirc/code/exe/$exef has been created."
echo "###########################"

