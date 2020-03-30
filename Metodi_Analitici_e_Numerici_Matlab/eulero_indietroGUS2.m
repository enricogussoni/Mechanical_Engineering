function [t_h,u_h,iter_nwt] = eulero_indietroGUS2(f,df,t_max,y0,h)
% ERRATO

t_h=0:h:t_max;
N=length(t_h);
u_h=zeros(1,N);
u_h(1)=y0;

% Parametri di Newton
tol=1e-12;
maxit=500;
iter_nwt=zeros(1,N); % Iterazioni di Newton ad ogni istante temporale

for i=2:N
    uv = u_h(i-1); % u all'istante t
    
    % Implementazione metodo di Newton
    fun=f(t_h);
    dfun=df(t_h);
    
    err=1;
    x=zeros(1,maxit);
    x(1)=uv;
    k=0;
    
    while (err>tol && k<maxit) 
        xv = x(i-1);
        xn = xv - fun/dfun; % si può aggiungere controllo di 'break' su df=0
        err=abs(xn-xv);
        x(i)=xn;
        k=k+1;
    end
    
    iter_nwt(i)=k;
    ft = x(end); % f(t,u) all'istante t+1
    un = uv + h*ft; % u all'istante t+1
    u(i)=un;
end
    
    