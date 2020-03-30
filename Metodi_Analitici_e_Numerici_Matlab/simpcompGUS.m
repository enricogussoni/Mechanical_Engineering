function I = simpcompGUS(a,b,N,fun)
% Implementa la formula di quadratura di Simpson composita su intervalli 
% equispaziati.

H=(b-a)/N;
nodi=a:H:b;
m=(a+H/2):H:(b-H/2);
I=(H/6)*sum([fun(nodi(1)) 2*sum(fun(nodi(2:N))) 4*sum(fun(m)) fun(nodi(end))]);