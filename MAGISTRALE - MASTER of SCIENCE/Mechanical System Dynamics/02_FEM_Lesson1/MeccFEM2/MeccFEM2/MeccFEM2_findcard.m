function MeccFEM2_findcard(fid,card)
%MECCFEM2_FINDCARD Look inside the specified file for a given string
%   MECCFEM2_FINDCARD(FID,CARD) scroll through the file identified by the
%   file-id FID for the string CARD and position the position indicator at
%   the beginning of the following line

maxiter=1e5;

frewind(fid);

while ~feof(fid)
   line=fgetl(fid);
   if ~isempty(strfind(line,card))
      return
   end
end

error('Error\nCard %s can not be found',card)

