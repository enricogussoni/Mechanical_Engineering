clc
close all
clear all

%% 1

f = @(x) atan(7*(x-pi/2)) + sin((x-pi/2).^3);
x = linspace(-1, 6, 1000);
y = f(x);
figure(1);
plot(x,y)
title('f(x)=atan(7*(x-pi/2)) + sin((x-pi/2)^3)');
xlabel('x');
ylabel('y');
grid on
y0 = zeros(1000,1);
hold on
plot(x,y0)

grid on

%% 2

clc

toll = 1e-10;
nmax = 1000;

alpha = pi/2;

df = @(x) 7./( 1 + 49 * ( x-pi/2 ).^2 ) + 3 * (x-pi/2).^2 .* cos( (x-pi/2).^3 );
plot(x,df(x),'r')
if(df(alpha) ~= 0)
    disp('Radice semplice')
else
    disp('Radice multipla')
end

x01 = 1.5;
[xvect_1,it_1]=newton(x01,nmax,toll,f,df,1);

x02 = 4;
[xvect_2,it_2]=newton(x02,nmax,toll,f,df,1);

stimap(xvect_1);