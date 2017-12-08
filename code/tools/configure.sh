#-----------------------
# configure script.
# This script should be run before compiling.
# It enables default configuration.
# 2012/04/27, JBG.
#-----------------------

nb_arg=$#
# Checks the number of arguments
if [ "$nb_arg" != "0" ]; then
  echo "Error: usage: sh configure.sh (no argument)"
  exit
fi


echo "########################"
echo "# configure script."
echo "# "
echo "# Step 1: some tests"
echo "# "
echo "  "


##################
# Step 1: tests


# makedepend will be used to create the makefile.
# Test on its presence:

mkd=$(which makedepend)

if [ "$mkd" = "" ];
 then
   echo "test1: failure. makedepend has not been found in PATH"
   mkd_present="F"
  else
   echo "test1: success. makedepend is present in PATH: $mkd"
   mkd_present="T"
fi


# netcdf option

dfn_cdf=$(cat ./optfile | sed "s%^define_netcdf%@@%g" | grep "@" | cut -d= -f2)

# dfn_cdf="T" is the user wants to use netcdf output
# dfn_cdf="F" is the user don't

if [ "$dfn_cdf" = "T" ];
 then
  
     # nc-config will be used to link the netcdf libraries.
     # Test on its presence:
  ncc=$(which nc-config)
  #echo $ncc
  if [ "$ncc" = "" ];
   then
     echo "test2: failure. nc-config has not been found in PATH"
     ncc_present="F"
    else
     echo "test2: success. nc-config appears in PATH : $ncc"
     ncc_prsent="T"
  fi

 else
  ncc="F"

fi 


##################

# Crossroads:
# If { makedepend is not present OR {nc_config is not present AND the user wants netcdf} }, stop
#   else,continue.


if [ "$ncc" = "" ] || [ "$mkd" = "" ];
 then
   echo "At least one test failed."
   echo "stop."
   exit
fi

echo "  "
echo "# step 1 successful."
echo "# "

echo "# Step 2: configure optfile and namelist"
echo ""



##################
# Step 2: 

#
# Configuration of optfile
#


echo "optfile is completed with"

opt_temp="_optfile.tmp"
opt_temp2="_optfile2.tmp"
cp ./optfile $opt_temp

dirc=$(pwd | sed "s%/code$%%g")

# The main directory (where the code files are) is written in $opt_temp:
echo "  - the main directory       : " $dirc
cat $opt_temp | sed "s/^dirc/@@/g" | cut -d@ -f1-2 | sed "s%@%dirc=$dirc%g" > $opt_temp2
mv $opt_temp2 $opt_temp


# Then, netcdf libraries paths:
if [ "$dfn_cdf" = "T" ];
 then
 
  netcdf_libs=$(nc-config --flibs)
  netcdf_incl=$(nc-config --includedir)
  echo "  - the netcdf library paths : " $netcdf_libs
  echo "  - the netcdf include paths : " $netcdf_incl
  echo "  - netcdf version           : " `nc-config --version`

  # libraries paths:
  netcdf_lib1=$(echo $netcdf_libs | sed "s/-l/@/g" | cut -d@ -f1 | sed "s/@/-l/g")
  cat $opt_temp | sed "s/^netcdf_dir_lib/@@/g" | cut -d@ -f1-2 | sed "s%@%netcdf_dir_lib=$netcdf_lib1%g" > $opt_temp2
  mv $opt_temp2 $opt_temp

  #libraries names:
  netcdf_lib2=$(echo $netcdf_libs | sed "s/-l/@@/g" | cut -d@ -f2- | sed "s/@@/-l/g;s/@/-l/g")
  cat $opt_temp | sed "s/^netcdf_lnk_flag/@@/g" | cut -d@ -f1-2 | sed "s%@%netcdf_lnk_flag=$netcdf_lib2%g" > $opt_temp2

  #include paths:
  cat $opt_temp2 | sed "s/^netcdf_dir_inc/@@/g" | cut -d@ -f1-2 | sed "s%@%netcdf_dir_inc=-I$netcdf_incl%g" > $opt_temp

fi

mv optfile _optfile.bak
mv $opt_temp optfile
rm -f $opt_temp2
rm -f _optfile.bak


#
# Configuration of namelist
#

# The only thing that is done is to direct the model output towards the default directory ../output.
echo " "

dirout=$(pwd | sed "s%code$%output/%g")
echo "namelist is completed with" 
echo "  - the default directory for the model output: $dirout"

cat namelist | sed "s%^dirout=%@@%g" | cut -d@ -f1-2 | sed "s%@%dirout=\"$dirout\",%g" > $opt_temp
mv namelist _namelist.bak; mv $opt_temp namelist; rm -f _namelist.bak

# End of step 2.
#######################


echo ""
echo "# end of configure script"
echo "########################"
