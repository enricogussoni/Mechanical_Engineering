function I = trapcompGUS(a,b,N,fun)
% Implementa la formula di quadratura del trapezio composita su intervalli 
% equispaziati.

H=(b-a)/N;
nodi=a:H:b;
I=H*sum([0.5*fun(nodi(1)) fun(nodi(2:N)) 0.5*fun(nodi(end))]);
