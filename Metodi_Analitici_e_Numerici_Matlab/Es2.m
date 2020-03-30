% Lab 12 Es2
% Dal tema d'esame del 18 febbraio 2015

% u" - u' = -x , x€[0,1],
% u'(0)=-1
% u(1)=1

% Problema misto Neumann-Dirichlet senza il termine di reazione
% Scegliendo un rilevamento costante questo non varierà la forma del
% problema (le sue derivate I e II saranno nulle)

% a=-1;
% b=-1;

% f = @(x) -x;

% In questo caso può veler la pena di cambiare segno all'equazione per 
% mettersi nella condizione standard:
% -u" + u' = x , x€[0,1]
% (Aiuta ad evitare errori nel calcolo dell'integrale p.p.)

a=1;
b=1;

f = @(x) x;
F = h*f(xdof);
F(1) = h/2 * f(x0) + 1;

...
