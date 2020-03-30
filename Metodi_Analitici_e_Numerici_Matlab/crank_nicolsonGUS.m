function [t_h,u_h,iter_nwt]=crank_nicolsonGUS(f,df,t_max,y0,h)
%ERRATO

t_h=0:h:t_max;
N=length(t_h);
u_h=zeros(1,N);
u_h(1)=y0;

maxit=100;
iter_nwt=zeros(1,maxit);
tol=1e-12;

for i=2:N
    uv = u_h(i-1);
    
    fun = @(t_h,u) u(i)+0.5*h*(f(t_h(i-1),u(i-1))+ f(t_h(i),u(i+1))) -u;
    dfun = @(t_h,u) 0.5*h*(df(t_h(i-1),u(i-1))+ df(t_h(i),u(i+1)))-1;
    
    err=1;
    k=0;
    xv = uv;
    
    while (k<maxit && err>tol)
           xn = xv + fun(xv)/dfun(xv);
           err=abs(xn-xv);
           k=k+1;
           xv=xn;
    end
    
    u_h(i)=xn;
    iter_nwt(i)=k;
    
end
           
        