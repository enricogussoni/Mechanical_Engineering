function [t_h,u_h] = eulero_avantiGUS(f,t_max,y0,h)

u=zeros(1,t_max);
u(1)=y0;
t_h=y0:h:t_max;
for i=2:length(t_h)
    u(i)=u(i-1)+h*f(t(i-1),u(i-1));
end