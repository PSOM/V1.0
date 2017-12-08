############################
# compile_complex.sh
#
# This script enables the combined use of model and model_build 
#  to create a hybrid executable with both old and new files.
# --------------------------
#
# V1.0 2012/02/27. JBG.
# V2.0 2013/01/25. JBG.
# V2.1 2013/03/06. JBG.  


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
echo "#"
}

######################
# Arguments

nb_arg=$#
rep_alt=$1


if [ "$rep_alt" = "-help" ] || [ "$rep_alt" = "--help" ]; then

cat >_compile.sh.help <<EOF

COMPILE.SH -- V2.1, 2013/03/06.
 
Usage: sh compile*.sh directory
 
Various cases of errors: 

- Error: usage: sh compile*.sh subdirectory.
     -> compile.sh is used in two different ways:
          - without argument. In that case, expe_template is used as a superceding directory. Because this directory is empty by default, the default model is compiled.
          - with an argument. In that case, the argument should be an existing subdirectory of ./code. This subdirectory will be used an "experiment" directory.
 
- Error: X is not a directory.
     -> the argument of compile must be an existing subdirectory of ./code.
     -> Solution: create a subdirectory following the procedure indicated in GET_STARTED_2
 
- Error: The last compilation was the experiment X.
  Error: You need to clean the environment before compiling the Y experiment.

     -> The user first compiled the experiment X and then tried to compile the experiment Y.
        The current structure of psom allows to develop several versions of the model. Nevertheless, to switch from one version to the other, you need to clean the directories and files needed to compile.
     -> Solution: run sh tools/clean.sh and re-try.

EOF

cat _compile.sh.help | sed "s/^/# /g"
rm -f _compile.sh.help
exit
fi



############################


echo "########################"
echo "# sh compile.sh"
echo "# "
echo "#----------------------------------------------------------------|"
echo "# INITIALIZATION                                                  "

# --- Checks -----

# Number of arguments
if [ "$nb_arg" != "0" ] && [ "$nb_arg" != "1" ]; then
  echo "Error: usage: sh compile*.sh subdirectory"
  exit
fi

if [ "$nb_arg" = "0" ]; then  # If no argument, expe_template is default.
  echo "# Default setting."
  dir_alt="expe_template"
  exe_ext=""
 else
  echo "# Custom setting."    # If an argument exists, it will be the alternative folder.
  dir_alt=$rep_alt
  exe_ext="_"$dir_alt
fi

echo "# "

# Check the nature of the candidate alternative folder: it has to be a directory and not be some of the key directories.
if [ "$dir_alt" == "doc" ] || [ "$dir_alt" == "examples" ] || [ "$dir_alt" == "exe" ] || [ "$dir_alt" == "mkfile" ] || [ "$dir_alt" == "model" ] || [ "$dir_alt" == "model_build" ]; then
  echo "Error: $dir_alt cannot be used as an experiment directory"
  exit
fi

if [ "$dir_alt" == "tools" ] || [ "$dir_alt" == "utils" ]; then
  echo "Error: $dir_alt cannot be used as an experiment directory"
  exit
fi

if [ ! -d "$dir_alt" ]; then
  echo "Error: $dir_alt is not a directory."
  exit
fi

######################
# Variables are initialized

optfile_gene="./optfile"                                              # optfile_gene is ./code/optfile
dirc=$(cat $optfile_gene | grep -v "#" | grep "dirc=" | cut -d= -f2)  # dirc is the ./code directory
dirloc=`pwd`
EXE="nh"
exef="$EXE""$exe_ext"                                                 # exef is the name of the executable that will be created

echo "# Use of model and $dir_alt to create $exef"
echo "# "

dirh=$dirc"/code/model"

# Default directories model/src and inc
dirsrc=$(echo $dirh"/src")
dirinc=$(echo $dirh"/inc")

# Alternative directories $dir_alt/src and inc
dirsrcalt=$(echo $dirh"/../$dir_alt/src")
dirincalt=$(echo $dirh"/../$dir_alt/inc")

# Compiling directories model/src and inc
dirsrcbui=$(echo $dirh"/../model_build/src")
dirincbui=$(echo $dirh"/../model_build/inc")

# Optfiles
optfile_add=$(echo $dirsrcalt"/../optfile_add")
optfile=$(echo $dirsrcbui"/../optfile")


######################
# Reading the potential model_build optfile

# This step gets info about the previous compilation attempt.

diralt_prev=""

# If the previous compilation ended correctly, continue, else, stop. 
if [ -f $optfile ]; then
  compstat_prev=$(cat $optfile | grep "compil_status" | cut -d= -f2)
  if [ "$compstat_prev" != "OK" ]; then
    echo "Error: Last compilation attempt ended prematurely."
    echo "Error: You have to clean the environment before trying to compile again."
    exit
  fi
fi


#If the previous compilation was the creation of another experiment, this case is not supported. 
# There is indeed a very high risk to mix files from one experiment to the other. 
if [ -f $optfile ]; then
  diralt_prev=$(cat $optfile | grep "dir_alt=" | cut -d= -f2)
  if [ "$diralt_prev" = "$dir_alt" ]; then  # The previous compilation dealt with the same experiment: it is safe to go forward with the same files in model_build.
    echo "# Compiling the same experiment..."
   else
    echo "Error: The last compilation was the experiment $diralt_prev ."
    echo "Error: You need to clean the environment before compiling the $dir_alt experiment."
    exit
  fi
 else
  echo "# Model_build is clean."
fi
#Note: we could systematically clean model_build to avoid this step. 



#######################
# Addition of objects


# The object list is the combination of the list in ./code/optfile and in ./$dir_alt/optfile_add .

# ./code/optfile
cat $optfile_gene | sed "s/^[ ]*#/@/g" | grep -v "@" > $optfile

# ./$dir_alt/optfile_add
if [ -f $optfile_add ]; then
  cat $optfile_add | sed "s/^[ ]*#/@/g" | grep -v "@" >> $optfile
fi

########################
# Addition of information about the current compilation

# the name of the experiment
echo " " >> $optfile
echo "dir_alt="$dir_alt >> $optfile

# the status of the current compilation
echo " " >> $optfile
echo "compil_status=IN PROGRESS" >> $optfile


# END OF THE HEADER
################################################
# START OF THE SCRIPT


##############
# Step1 : Detecting alternative files.

# List the alternative directories
lsrcalt=$(ls $dirsrcalt)
lincalt=$(ls $dirincalt)

list_of_files

echo "# +------                                                         "
echo "# |                                                               "
echo "# | The following files are detected in the alternative directory:"
echo "# |  /src : " $lsrcalt
echo "# |  /inc : " $lincalt
echo "# |                                                               "
echo "# +------                                                         "


##############
# Step2: Copying the files.

# Any file in an alternative directory, even a very old one, will be considered as recent.

# Copying the files from "model" with their status (including date)
cp -fp $dirsrc/* $dirsrcbui 2> cmp.tmp
cp -fp $dirinc/* $dirincbui 2>> cmp.tmp


# In order to match the "makefile" philosophy (based on dates), every file found in the alternative directories have their timestamp updated to the current time.
# Copy the files from the alternative directory with the current date.

if [ "$lsrcalt" != "" ]; then
 cp $dirsrcalt/* $dirsrcbui
fi


if [ "$lincalt" != "" ]; then
cp $dirincalt/* $dirincbui
fi

# As a security, the copied files are touched.
for fic in $lsrcalt
 do
  touch $dirsrcbui/$fic
done

for fic in $lincalt
 do
  touch $dirincbui/$fic
done


# For every file that is found in the alternative src directory, we check its presence in the list of objects in ./optfile.
# It is only informative.

  bobjf=""
  bobjn=""

for srcalt in $lsrcalt
 do
  srcalto=$(echo "@"$srcalt"@" | sed "s/\.f90//g")
  srcalto2=$(echo $srcalto | sed "s/@//g")
  sloc=$(cat $optfile | grep "="  | sed "s/=/=@/g;s/$/@/g" | sed "s/ /@/g" | sed "s/lobj/X@@X/g;s/lmod/X@@X/g" | sed "s/\.mod//g;s/\.o//g" | grep "X@@X" | cut -d= -f2 )
  sloc=$(echo $sloc | grep $srcalto)

  if [ "$sloc" = "" ]; then
    bobjn=$(echo $bobjn $srcalt)
    # echo "Warning ! " $srcalto2 " does not appears in the list of objects nor in the list of modules in optfile. $srcalto2 will not be taken into account during compilation !!"
   else
    bobjf=$(echo $bobjf $srcalt)
    #echo $srcalto2 "appears in the list of objects or modules in optfile. It will be taken into account during compilation."
  fi

done 


nobj=$(cat $optfile | grep "="  | sed "s/=/=@/g;s/$/@/g" | sed "s/ /@/g" | sed "s/lobj/X@@X/g;s/lmod/X@@X/g" | sed "s/\.mod//g;s/\.o//g" | grep "X@@X" | cut -d= -f2  | sed "s/@/ /g" | wc -w)

echo "# "
echo "# The list of objects contains $nobj objects. " | sed "s/ [ ]*/ /g" 

 
if [ "$bobjf" != "" ]; then
 echo "#  "
 echo "# +------                                                         "
 echo "# |"
 echo "# | The following files appear in the list of objects or modules in optfile. They will be taken into account during compilation: "
 echo "# |" 
 echo "# |     " $bobjf
 echo "# |" 
 echo "# +------                                                         "
fi

if [ "$bobjn" != "" ]; then
 echo "# "
 echo "# +------                                                         "
 echo "# |"
 echo "# | The following files do NOT appear in the list of objects or modules in optfile. They will NOT be taken into account during compilation:"
 echo "# |"
 echo "# |     " $bobjn
 echo "# |"
 echo "# +------                                                         "
fi

list_of_files


##############
# Step3 : Compiling

echo "# "
echo "# END OF THE INITIALIZATION"
echo "#----------------------------------------------------------------|"
echo "# CREATION OF THE MAKEFILE "

# Calling genmakefilel creates the makefile based on the information stored in model_build/optfile
$dirc/code/tools/genmakefilel $optfile > cmp.tmp

cat cmp.tmp | sed "s/^/# /g" > cmp2.tmp; cat cmp2.tmp  # output

echo "# "
echo "#----------------------------------------------------------------|"
echo "# BUILDING THE EXECUTABLE"
echo "# "

# Building of the executable
make -C $dirc/code/mkfile $EXE
#make -C $dirc/code/mkfile clean
cp $dirc/code/mkfile/$EXE $dirc/code/exe/$exef

echo "#-------------------------------------------------------"


###############
# Step4: Ending 

list_of_files

echo "# "
echo "# `ls -al $dirc/code/exe/$exef`"

rm -f cmp.tmp cmp2.tmp

# model_build/optfile is updated with the new status of the compilation ("IN PROGRESS" -> "OK")
cat $optfile | sed "s/compil_status=IN PROGRESS/compil_status=OK/g" > cmp.tmp; mv cmp.tmp $optfile; 

echo "###########################"

