x0=0;
xl=1;
h=0.1;
xdof=[x0:h:xl-h]'; 
N=length(xdof);
U=ones(N,1);
U1=ones(N-1,1);
a=1; b=-1;
A=a/h*(2*diag(U)-diag(U1,1)-diag(U1,-1))+b*(0*diag(U)-0.5*diag(U1,-1)+0.5*diag(U1,1));
A(1,1)=a/h+b/2;   % modifica del termine A(1,1) per dato di Neumann in x=0.

frhs=@(x)-x;
f=h*frhs(xdof);
f(1)=h/2*frhs(x0)+1; % modifica del termine f(1) per dato di Neumann in x=0. Attenzione ai segni!!

u_tilde=A\f; 
xbc=[xdof;xl]; 
ubc_tilde=[u_tilde;0]; 
ubc = ubc_tilde+1;  %alla soluzione del problema con dato omogeneo, sommo la funzione rilevamento, 
%ovvero la funzione costante g(x)=1
plot(xbc,ubc);
legend('soluzione elementi finiti');