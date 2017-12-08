function varanom=var_anom(iw,var)

varmean=runmean(var,iw,3);
varanom=var-varmean;

