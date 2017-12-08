# Makefile must have been created and in code/mkfile/.

EXE="nh"

sh tools/genmakefilel
make -Cmkfile/ $EXE
cp ./mkfile/$EXE ./exe/$EXE

