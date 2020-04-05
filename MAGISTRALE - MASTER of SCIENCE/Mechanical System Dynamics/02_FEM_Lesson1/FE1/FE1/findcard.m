function tmp=findcard(fid,card)
% Cerca all'interno del file specificato da fid la stringa 'stringa' e posiziona il puntatore di lettura file alla riga successiva.
% Restituisce 1 se la stringa e' stata trovata, 0 altrimenti.
maxiter=1e5;

frewind(fid);
for i=1:maxiter
   if feof(fid)
      error(['Impossibile trovare stringa ' card])
   end
   riga=fgets(fid);
   if ~isempty(findstr(riga,card))
      return
   end
end

error(['Impossibile trovare stringa ' card])

