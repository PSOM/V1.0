function gradup=presswork(anom)

  global dir file nt iy0 iyf iz0 izf
  global anom ntanom
  global www vvv yyy zzz

  file=[dir 'xslice_005.cdf'];

  Jac=ncread(file,'Jac',[iy0 iz0], [iyf-iy0 izf-iz0]);

  file=[dir 'xslice_face_005.cdf'];

  vflux=ncread(file,'vf'  ,[iy0 iz0 1],[iyf-iy0+1 izf-iz0 nt]);
  wflux=ncread(file,'wf'  ,[iy0 iz0 1],[iyf-iy0   izf-iz0 nt]);
  ptfy =ncread(file,'ptfy',[iy0 iz0 1],[iyf-iy0+1 izf-iz0 nt]);
  ptfz =ncread(file,'ptfz',[iy0 iz0 1],[iyf-iy0   izf-iz0 nt]);

  nyl=size(vflux,1);
  nzl=size(wflux,2);

  if anom==1
    for it=1:nt
      vwind (:,:,it)=ones(nyl  ,1)*vflux(1,:,it);
      wwind (:,:,it)=ones(nyl-1,1)*wflux(1,:,it);
      pywind(:,:,it)=ones(nyl  ,1)*ptfy (1,:,it);
      pzwind(:,:,it)=ones(nyl-1,1)*ptfz (1,:,it);
    end
    vfluxdat=var_anom(ntanom,vflux);
    wfluxdat=var_anom(ntanom,wflux);
    pfydat  =var_anom(ntanom,ptfy) ;
    pfzdat  =var_anom(ntanom,ptfz) ;
  elseif anom==2
    for it=1:nt
      dum  =squeeze(vflux(1,:,it));
      dumM =ones(nyl,1)*dum;
      vfluxdat(:,:,it)=vflux(:,:,it)-dumM;
      clear dumM dum

      dum  =squeeze(wflux(1,:,it));
      dumM =ones(nyl-1,1)*dum;
      wfluxdat(:,:,it)=wflux(:,:,it)-dumM;
      clear dumM

      dum  =squeeze(ptfy(1,:,it));
      dumM =ones(nyl,1)*dum;
      pfydat(:,:,it)=ptfy(:,:,it)-dumM-ptfy(:,:,1);
      clear dumM

      dum  =squeeze(ptfz(1,:,it));
      dumM =ones(nyl-1,1)*dum;
      pfzdat(:,:,it)=ptfz(:,:,it)-dumM-ptfz(:,:,1);
      clear dumM
    end
  else
    vfluxdat=vflux;
    wfluxdat=wflux;
    pfydat  =var_anom(ntanom,ptfy) ;
    pfzdat  =var_anom(ntanom,ptfz) ;
  end

  for it=1:nt
  for j =1:nyl-1
  for k =1:nzl-1
    dvfpdy=vfluxdat(j+1,  k,it)*pfydat(j+1,  k,it)-vfluxdat(j,k,it)*pfydat(j,k,it);
    dwfpdz=wfluxdat(  j,k+1,it)*pfzdat(  j,k+1,it)-wfluxdat(j,k,it)*pfzdat(j,k,it);

    gradup(j,k,it)=-(dvfpdy+dwfpdz)/Jac(j,k)/1.d5;
  end
  end
  end

end
