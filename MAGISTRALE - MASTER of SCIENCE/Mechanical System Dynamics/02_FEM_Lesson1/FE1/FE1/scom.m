function line=scom(fid)
line=fgetl(fid);
tmp=sscanf(line,'%s',1);
while(tmp=='!')
    if feof(fid)
      error('Impossibile individuare riga non commentata')
   end
 line=fgetl(fid);
 tmp=sscanf(line,'%s',1);
%  disp('in scom')
%  pause
end 

