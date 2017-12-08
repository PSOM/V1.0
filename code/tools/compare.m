cd /Applications/
addpath ([pwd '/mexcdf/mexnc']); 
addpath ([pwd '/mexcdf/snctools']); 
javaaddpath([pwd '/mexcdf/netcdfAll-4.2.jar']);

set(gcf,'PaperUnits','centimeters')
xSize = 22; ySize = 26;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[0 0 xSize*50 ySize*50])

get(0,'Screensize');set(0,'Units','normalized');% set pop up figure
                                                % page size

                                                
                                               
%%%%%%%%%%%
% go to the folder in which this script is.
cd /Users/jgilet/work_on_structure/V0.5_2nd_trial/code

% choose the number of time steps saved in ktopsurfmovie files
%nt=1001;

nt=6;
%nt=size(double(nc_varget(filename,'u')))



% choose the option:
% option=x : ktopsurf stored in kth choice 
% option2=1 : just visualize the files
% option2=2 : make differences between these files, with VARG.

%option=1;
option2=2;

option_sauv=0                                     
                                                
%%%%%%%%%%%

if option2==1
  if nt==1001
    loption=1:3  %number of files that the script reads.
    ltimesteps=1:100:1001 %series of time steps that are visualized.
   else
    % loption=1:2
    loption=[1,2]
    ltimesteps=nt %:10:nt    
  end
 else
  loption=1
  ltimesteps=1:nt %10:nt %10:nt %10:nt %10:nt %2 %nt:nt % 10:nt %1:10:101 
end

for option=loption
       
  % selection of the file that is read
    switch option
  %    case 1;  filename2='zslice_016.cdf'; filename1='/Users/jgilet/PDC/V0.4/output_without_options/'
  %    case 2;  filename2='zslice_016.cdf'; filename1='/Users/jgilet/PDC/V0.4/output/'

   %   case 1;  filename2='gslice_005.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_106B_gfortran/'
   %   case 2;  filename2='gslice_005.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_106B_pgf95/'
   
   %   case 1;  filename2='zslice_032.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_wiggle/'
   %   case 2;  filename2='zslice_032.cdf'; filename1='/Users/jgilet/work_on_structure/V0.41/output/'
   
   %   case 1;  filename2='zslice_032.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_wiggle/'
   %   case 2;  filename2='zslice_032.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_wiggleT/'
    
   %   case 1;  filename2='zslice_032.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_wiggle_FL_0/'
   %   case 2;  filename2='zslice_032.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_wiggle_FL_1/'
    
          
   %   case 1;  filename2='gslice_090.cdf'; filename1='/Users/jgilet/work_on_structure/V0.4_Sche/output_106B_0p_1p_1p_1p_2000p_11p_1p_1p/'
   %   case 2;  filename2='gslice_090.cdf'; filename1='/Users/jgilet/work_on_structure/V0.4_Sche/output_106B_0p_1p_1p_1p_2000p_11p_1p_1p/'

   
   %   case 1;  filename2='zslice_032.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_wiggle/'
   %   case 2;  filename2='zslice_032.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_wiggle_FL/'

   
    %  case 1;  filename2='zslice_032.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_wiggle_FL/'
    %  case 2;  filename2='zslice_032.cdf'; filename1='/Users/jgilet/work_on_structure/V0.41/output_3500/'
   
    %  case 1;  filename2='zslice_032.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_NA_test/'
    %  case 2;  filename2='zslice_032.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/code/melissa_NA/output/'  
    
    %  case 1;  filename2='zslice_001.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_wiggle_FL_2/'
    %  case 2;  filename2='zslice_001.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_wiggle_FL/'

    %  case 1;  filename2='zslice_001.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_106B/'
    %  case 2;  filename2='zslice_001.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_106B_2/'

    %  case 1;  filename2='zslice_032.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_106B/'
    %  case 2;  filename2='zslice_032.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_test_advect/'

    
    %  case 1;  filename2='zslice_016.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_wiggle_FL_pickup/PICKUP_0/'
    %  case 2;  filename2='zslice_016.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_wiggle_FL_pickup/PICKUP_1/'

     case 1;  filename2='gslice_005.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_wiggle_FL_10_un_jet/'
     case 2;  filename2='gslice_005.cdf'; filename1='/Users/jgilet/work_on_structure/V0.5_2nd_trial/output_wiggle_FL/'

    
    end
    

  filename=strcat(filename1,filename2);

  if option2==1                 
    % Reading the file                                      
    %xc=double(nc_varget(filename,'xc')); yc=double(nc_varget(filename,'yc'));
    xc=double(nc_varget(filename,'xc')); yct=double(nc_varget(filename,'yc'));yc=yct(:,1)    
    nx=length(xc);ny=length(yc);
    VAR=zeros(nt,ny,nx,8);
    VAR(:,:,:,1)=double(nc_varget(filename,'u'));
    VAR(:,:,:,2)=double(nc_varget(filename,'v'));
    VAR(:,:,:,3)=double(nc_varget(filename,'w'));
    VAR(:,:,:,4)=double(nc_varget(filename,'rho'))-1000.;
  end

  
  % We prepare for a possible future call of option2=2
  if option2==1
    %VARG=zeros(nt,ny,nx,8,4);  % Sloppy coding...
    for k=1:4
      VARG(:,:,:,k,option)=VAR(:,:,:,k);
    end
  end

  % definition of the names of the variables
  nvar=['u','v','w','rho'];


  for l=ltimesteps
    
    for k=1:4
    
      grid on

      %initanim

      %----- examples ---------
      %text(-50,180,'wind')
      %axis off
      %cl=colorbar('EastOutside');
      %set(cl,'YAxisLocation','right','FontSize', 11)
       %set(ps,'Color',[1 1 1]) 
      %set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);

      %----- surface view------

      % horizontal box visualization   
      %switch k
      %  case 1
      %    xmin=0.12;ymin=0.76;dx=0.66;dy=0.23;xmin2=0.91;ymin2=ymin;dx2=0.04;dy2=dy;
      % case 2
      %    xmin=0.12;ymin=0.51;dx=0.66;dy=0.23;xmin2=0.91;ymin2=ymin;dx2=0.04;dy2=dy;
      %  case 3
      %    xmin=0.12;ymin=0.26;dx=0.66;dy=0.23;xmin2=0.91;ymin2=ymin;dx2=0.04;dy2=dy;
      %  otherwise
      %    xmin=0.12;ymin=0.01;dx=0.66;dy=0.23;xmin2=0.91;ymin2=ymin;dx2=0.04;dy2=dy;
      %end

  
      switch k
        case 1
          xmin=0.06;ymin=0.52;dx=0.36;dy=0.45;xmin2=0.45;ymin2=ymin;dx2=0.04;dy2=dy;
        case 2
          xmin=0.56;ymin=0.52;dx=0.36;dy=0.45;xmin2=0.94;ymin2=ymin;dx2=0.04;dy2=dy;
        case 3
          xmin=0.06;ymin=0.02;dx=0.36;dy=0.45;xmin2=0.45;ymin2=ymin;dx2=0.04;dy2=dy;
        otherwise
          xmin=0.56;ymin=0.02;dx=0.36;dy=0.45;xmin2=0.94;ymin2=ymin;dx2=0.04;dy2=dy;
      end
    
      if option2==1
        arplot=VAR(l,:,:,k);
      else
        arplot=VARG(l,:,:,k,2)-VARG(l,:,:,k,1);
      end

      if (option2==1) 
        if k<1.5
           varmin=-1.; varmax=1.; varinc=0.002;
         elseif k<2.5
           varmin=-0.1; varmax=0.1; varinc=0.002;
         elseif k<3.5
             varmin=-0.0001; varmax=0.0001; varinc=0.000005;
         else     
          varmoy=0.5*(max(max(arplot))+min(min(arplot)));vardel=max(0.5*(max(max(arplot))-min(min(arplot))),1e-13);
          varmin=varmoy-vardel;varmax=varmoy+vardel;
          varinc=0.02*(varmax-varmin);   
          % varmin=-0.5*(max(max(arplot))-min(min(arplot)));varmax=-varmin;varinc=0.02*(varmax-varmin);
        end
       else
            
        if k<3.5   
          varmin=-0.5*(max(max(arplot))-min(min(arplot)));varmax=-varmin;varinc=0.02*(varmax-varmin);
         else
          varmoy=0.5*(max(max(arplot))+min(min(arplot)));vardel=max(0.5*(max(max(arplot))-min(min(arplot))),1e-13);
          varmin=varmoy-vardel;varmax=varmoy+vardel;
          varinc=0.02*(varmax-varmin);   
        end
     end
  
      varint=varmin:varinc:varmax;



      % Plotting of the values
      h=figure(option+100*option2);
      subplot('Position',[xmin ymin dx dy]);
      pcolor(squeeze(xc),squeeze(yc),squeeze(arplot)); shading interp; caxis([varmin varmax]);
      %xlabel('dist (km)'): 
      ylabel('dist (km)');
      xt1=strcat(nvar(k),'    max : ',num2str(max(max(arplot))),',  min : ',num2str(min(min(arplot))))
      % In the title of the 3rd subplot, the timestep in included
      if k==3
        title(strcat(xt1,'     t=',num2str(l)))   
        % title(strcat(xt1,filename))
       else
        title([xt1])
      end
    
      ps= text(nx/10,ny/12,strcat(filename1));ps= text(nx/10,ny/20,strcat(filename2));
      subplot('Position',[xmin2 ymin2 dx2 dy2]);
      axis off;
      cl=colorbar('West');caxis([varmin varmax]);


    end % loop on k

    %outname=strcat('output_matlab_',num2str(l),'.eps')
    %print(h,'-dps',outname)
    if option_sauv==1
      % The plot is saved in png files.
      % output_matlab__(file)__(nb of time steps in ktop...)__(time step)
      if option2==1
        outname=strcat('output_matlab_',num2str(option),'__',num2str(nt),'__',num2str(l),'.png')
        print(h,'-dpng',outname)
      else
        outname=strcat('output_matlab_diff_5_4','__',num2str(nt),'__',num2str(l),'.png')       
        print(h,'-dpng',outname)        
      end    
      %arplot=reshape(abs(fft(VAR(l,ny/2,:,3))),[nx,1])';
      %figure(option+10);
      %  plot(1:nx/2,arplot(1:(nx/2)));
    end


    %MATLAB increase label size
    hAll = findall(gcf);
    for idx = 1 : length(hAll)
      try
        set(hAll(idx),'fontsize',10);
      catch
        % never mind...
      end
    end
    
    set(gcf,'defaultaxesfontsize',10)

    % saveframe(l)

  end % loop on l (time steps)
end % loop on option (files)
