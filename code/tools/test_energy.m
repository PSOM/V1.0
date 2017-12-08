

fid = fopen('tools/_values_reference'); Ar = fscanf(fid, '%g %g',[2 10]);fclose(fid);


fid = fopen('_values_current'); Ac = fscanf(fid, '%g %g',[2 10]);fclose(fid);

diff_max=max(max(abs(Ar-Ac)./(Ar+Ac)));


fid = fopen('_value_diff','w');fprintf(fid, '%g', diff_max);fclose(fid);
