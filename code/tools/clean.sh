# Makefile must have been created and in code/mkfile.


dt=$(date +"%m_%d_%Y_%H_%m_%S")
[ -d model_build ] && mv model_build trash/model_build_$dt

mkdir -p model_build/src/ model_build/inc/

echo "clean model_build"
echo "clean mkfile"

make -Cmkfile/ clean


