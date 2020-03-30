% Esercizio strutture dati

% Puliamo il workspace
clear 
clc

arrayProfessori = struct;
arrayProfessori(1).matricola = 123123;
arrayProfessori(1).nome = 'Wanna ';
arrayProfessori(1).secondoNome = '';
arrayProfessori(1).cognome = 'Marchi';
arrayProfessori(1).corsi = struct; % ci limitiamo a dire che e' una struct
whos
disp(arrayProfessori) % visualizziamo la struttura di arrayStudenti
                      % non il contenuto. 

arrayProfessori(2).nome = 'Andrea';
arrayProfessori(2).secondoNome % attenzione! e' una lettura!
                               % il campo c'e' perche' un'altro elemento 
                               % dell'array di strutture, specificamente
                               % arrayProfessori(1), ha definito quel campo
                               % che per questo motivo lo ha ereditato
                               % Riportato come ans = [], size 0x0
disp('Riportato come ans = [], size 0x0; inizializzato perche qualche')
disp('altra struttura aveva questo campo, ma qui non e stata assegnata')
arrayProfessori(2).cognome = 'Dipre';
arrayProfessori(2).corsi = struct;
whos
disp('Perche non il doppio di 912? Perche quando istanzio una struttura')
disp('dati la memoria effettivamente usata e quella dei dati piu un po')
disp('di memoria per gestire la struttura stessa; si paga la prima volta,')
disp(' poi cresce piu o meno come i dati che salviamo dentro')

disp(arrayProfessori) % visualizziamo la struttura di arrayProfessori
                      % non il contenuto. 
var = arrayProfessori(2).corsi;
whos
whos var

% ora inseriamo i corsi
% aggiungiamo DINAMICAMENTE due campi, nome e urlCorso
% osservazre come una struttura possa diventare un array di strutture
% semplicemente spiazzando la variabile con un indice
arrayProfessori(2).corsi(1).nome = 'La Vita e le Opere del Maestro Osvaldo Paniccia';
arrayProfessori(2).corsi(1).urlCorso = 'goo.gl/Z5DlU4';
arrayProfessori(1).corsi(1).nome = 'La Chimica del Sale del Mar Morto';
arrayProfessori(1).corsi(1).urlCorso = 'goo.gl/a95qrU';
arrayProfessori(1).corsi(2).nome = 'Estetica della Comunicazione 1';
arrayProfessori(1).corsi(2).urlCorso = 'goo.gl/EOIAWQ';
whos
disp(arrayProfessori) % visualizziamo la struttura di arrayProfessori
                      % non il contenuto. 

% Possiamo anche aggiungere dinamicamente un campo struttura (creata dyn)
arrayProfessori(2).corsi(1).distribuzioneOre = struct;
arrayProfessori(2).corsi(1).distribuzioneOre.lezioneFrontale = 100;
arrayProfessori(2).corsi(1).distribuzioneOre.esercitazione = 10;
arrayProfessori(2).corsi(1).distribuzioneOre.laboratorio = 1000;
whos
disp(arrayProfessori) % visualizziamo la struttura di arrayProfessori
                      % non il contenuto. 

disp ('I corsi del primo professore (non possiamo usare disp perche e un array)')
arrayProfessori(1).corsi.nome
disp ('Il primo corso del primo professore (tutti i campi)')
disp(arrayProfessori(1).corsi(1))
disp ('I corsi di tutti i professori')

% arrayProfessori.corsi.nome? 
% NO
% Errore: Dot name reference on non-scalar structure. 
% Error in studenti (line 68)
% arrayProfessori.corsi.nome

% arrayProfessori(:).corsi.nome?
% NO
% Scalar index required for this type of multi-level indexing.
% Error in studenti (line 75)
% arrayProfessori(:).corsi.nome

% arrayProfessori(:).corsi(:).nome?
% NO
% Scalar index required for this type of multi-level indexing.
% Error in studenti (line 81)
% arrayProfessori(:).corsi(:).nome

for i = 1:length(arrayProfessori)
    arrayProfessori(i).corsi.nome
end

% ATTENZIONE
% cosa restituisce disp (arrayProfessori.corsi.distribuzioneOre)?
% un errore: inatti quella ? "l'esplosione" di tutti i professori, tutti i
% corsi e tutte le distribuzioni orarie! Cosa che la disp non puo'
% visulaizzare
              
% ATTENZIONE
% leggere un campo prima che ne sia stata segnalata l'esistenza a MATLAB
% e' un ERRORE. PRIMA dovete "appiccicarlo" alla struttura dati, poi lo 
% potete leggere ANCHE SE NON LO AVETE INIZIALIZZATO, come alla riga 18.

% FUNZIONI AVANZATE MATLAB
% Se volete vedere i video direttamente da MATLAB, usate il comando
% url = ... % lo dovete recuperare dalla struct; come si fa?
% web(url, '-browser')