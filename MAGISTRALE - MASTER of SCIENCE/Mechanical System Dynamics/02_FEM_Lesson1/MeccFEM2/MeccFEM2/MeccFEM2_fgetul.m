function line=MeccFEM2_fgetul(fid)
%MECCFEM2_FGETUL Get the next un-commented line
%   LINE=MECCFEM2_FGETUL(FID) scroll through the file identified by the
%   file-id FID and returns the next line that does not begins with the
%   comment mark !



%Get the first line
line=fgetl(fid);

%Look at the first char
marker=sscanf(line,'%s',1);

%If it's the comment maker OR the line is empty -> get the next line
while isempty(marker) || marker(1)=='!'
    
    % If the eof is reached throw back an error
    if feof(fid)
      error('Reached the End of File')
    end
 
 line=fgetl(fid);
 marker=sscanf(line,'%s',1);
end 

