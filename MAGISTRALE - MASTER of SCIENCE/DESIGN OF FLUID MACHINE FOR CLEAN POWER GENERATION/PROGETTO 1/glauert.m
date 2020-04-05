function[a a_p] = glaubert(x) %x local speed ratio vettore
a0 = (0.2600:0.0001:0.3333)';
a_p0 = (1-3*a0)./(4*a0-1);
x0 = sqrt((a0.*(1-a0))./(a_p0.*(1+a_p0)));
a=interp1(x0,a0,x);
a_p= (1-3*a)./(4*a-1);
% figure()
%     plot(x,a,'o'); hold on;
%     plot (x0,a0,'o'); hold on;
% legend('a','a0');

end






