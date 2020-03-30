%es1
a=1;
b=0;
c=1;

x0=0;
xL=1;

u0=0;
u1=0;

%% 1
h=0.1;
xnodi=[x0+h:h:xL-h]'; % non servono i nodi estremi per la cond. di D.
n=length(xnodi);

Dl=ones(1,n);
Dc=ones(1,n-1);

Aa = 2*diag(Dl) - diag(Dc,1) - diag(Dc,-1);
Ac = 2/3 *diag(Dl) + 1/6 * diag(Dc,1) + 1/6 * diag(Dc,-1);
A = a/h * Aa + c*h * Ac;

%%
F = h*(-1) + 0.*xnodi;

%% 2
if (min(eig(A))<0)
    fprintf('A non definita positiva\n');
else
    fprintf('A definita positiva\n');
end

%% 3
k=cond(A);

%% 4
u=pcg(A,F);

%% 5
u_es = @(x) -1-((2*cosh(1)-1)*sinh(x)-sinh(x+1))/sinh(1);
xdis=linspace(x0,xL,1000);

xdef=[0;xnodi;xL];
udef=[0;u;0];

plot(xdis,u_es(xdis),xnodi,u)

%% 6
integranda = @(x) (u_es(x) - interp1(xdef,udef,x)).^2;
normL2 = sqrt(simpcomp(x0,xL,1000,integranda));
