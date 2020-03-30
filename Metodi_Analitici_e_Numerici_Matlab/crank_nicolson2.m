function[th,uh,itN]=crank_nicolson2(f,df,T,y0,h)

th=0:h:T;
n=length(th);
uh=zeros(1,n);
uh(1)=y0;

tol=1e-6;
maxit=200;
itN=zeros(1,n);

for c=2:n
    uv=uh(c-1);
    tn=th(c);
    tv=th(c-1);
    
    fun = @(u) uv + h/2 * (f(tv,uv)+f(tn,u))-u;
    dfun = @(u) h/2 * df(tn,u)-1;
    
    err=1;
    it=0;
    xv=uv;
    while(err>tol && it<maxit)
        %if controllo derivata...
        xn = xv - fun(xv)/dfun(xv);
        err=abs(xn-xv);
        it=it+1;
        xv=xn;
    end
    uh(c)=xn;
    itN(c)=it;
end