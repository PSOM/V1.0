clear

% This code converts the binary file ouput from the model into CSV files.
% This can be useful to open in other programs (e.g., import in a SQL
% database)

%%
%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%

% number of files to import into matlab. If particles output split amongst
% several bin-files, this routine can loop through the incrementally named
% bin-files (incremental filenames is implemented by default in PSOM for
% bin-files)
number_of_files = 1;
% filepath of files to import - %'specify the path of the files to import here';
pathin = ~;
% Number of variables recorded (see particles.f90)
% WARNING !!!! %
% This number should be modified according to the parti_save
% subroutine in particles.f90. If the number of ouputs printed to the bin
% files is modified, the number below should be modified accordingly
number_of_variables = 20;
% filepath of the csv-file to be written - %'specify the path where you would like to write the csv-files';
pathout = ~;

%%%%%%%%%%%%%%%%%%%%%%%%
% CORE CODE
%%%%%%%%%%%%%%%%%%%%%%%%

% Loops through the number of files to open
for filenum = 1:number_of_files
    
    % Create fullpath
    filename = ['op.parti-',num2str(filenum,'%03.f'),'.bin'];
    fullpath = [pathin,'/',filename];
    
    % Display the file being extracted
    disp(['Converting ',filename,' into a CSV file ...']);
    
    % Open the file
    fileID = fopen(fullpath);
    
    % Extract the data (refer to particles.f90 to confirm that number)
    A = fread(fileID,[number_of_variables Inf],'double');
    A = A';
    
    % Remove all the zeros recorded
    ind = find(A(:,1)==0);
    if isempty(ind)~=1
        warning(['Missing ',num2str(length(ind)),' records...!'])
        A(ind,:) = [];
    end; clear ind
        
    % Write the output as a CSV file
    % Open the file
    fileID = fopen([pathout,'/newparticules-',num2str(filenum,'%03.f'),'.csv'],'w');
    % Write the data
    fprintf(fileID,'%d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',A');
    % Close the file
    fclose(fileID);
    
    % Clean-up
    clear A partnum fileID filename path fullpath
    
end; clear filenum