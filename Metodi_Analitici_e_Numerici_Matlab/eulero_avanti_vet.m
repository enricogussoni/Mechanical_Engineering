function [t_h,u_h]=eulero_avanti_vet(f,t_max,y_0,h)

t0=0;
t_h=t0:h:t_max;

N=length(t_h);
u_h=zeros(2,N);

u_h(:,1)=y_0;
for it=2:N
    u_old=u_h(:,it-1);
    u_h(:,it)=u_old+h*f(:,t_h(it-1),u_old);
end