
make -Cmkfile/ clean
make -Cmkfile/ nh 2> _log.tmp


file="_log.tmp"

nb_err=$(cat $file | grep "warning" | wc -l)


cat $file | grep "warning" | sed "s/warning/@/g" | cut -d@ -f2- | sort | uniq

echo "------------"
echo $nb_err
