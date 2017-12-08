# Makefile must have been created and in code/mkfile.


dt=$(date +"%m_%d_%Y_%H_%m_%S")
mv model_build trash/model_build_$dt

mkdir model_build/
mkdir model_build/src/
mkdir model_build/inc/

echo "clean model_build"
echo "clean mkfile"

make -Cmkfile/ clean


