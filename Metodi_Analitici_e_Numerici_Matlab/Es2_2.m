% Lab 9 Es 2
% L'equazione del calore

%1
a=1;
x=linspace(0,pi,1000);
f = @(x) sin(x) + 0.5*sin(2*x) + 0.2*sin(5*x);
N=100;
n=1:N;
t=[0,0.5,1];

xf=f(pi*n./(N+1));
yf=dst(xf);
bf=yf.*2/(N+1);

% u=zeros(length(t),length(n));
% for c=n
%     u(1)=sum(bf*exp(-a*t(1)*c^2).*sin(c*x));
%     u(2)=sum(bf*exp(-a*t(2)*c^2).*sin(c*x));
%     u(3)=sum(bf*exp(-a*t(3)*c^2).*sin(c*x));
% end

u = zeros(length(t),length(x));
for c=n
    u(1,:) = u(1,:) + bf(c)*exp(-a*t(1)*c^2).*sin(c*x);
    u(2,:) = u(2,:) + bf(c)*exp(-a*t(2)*c^2).*sin(c*x);
    u(3,:) = u(3,:) + bf(c)*exp(-a*t(3)*c^2).*sin(c*x);
end
 
figure(1)
plot(x,u(1,:),x,u(2,:),x,u(3,:))
legend('t=0','t=0.5','t=1')
title('Eq. del calore')
xlabel('coordinata spaziale')
ylabel('temperatura')

%2
t=0;
T=5;
h=0.05;
tempo=t:h:T;

figure(2)
for p=tempo
    u=0;
    for c=n
    u = u + bf(c)*exp(-a*p*c^2).*sin(c*x);
    end
    plot(x,u)
    hold on
    drawnow
    xlabel('coordinata spaziale')
    ylabel('temperatura')
    title('Andamento di T nel tempo')
end

%3
% f(x) = 1 if pi/4 < x < 3*pi/4 else 0

f = @(x) 1*(x>(pi/4)).*(x<(3/4)*pi);
T=1;
t=linspace(0,T,100);

xf=f(pi*n./(N+1));
yf=dst(xf);
bf=yf.*2/(N+1);

figure(3)
for p=tempo
    u=0;
    for c=n
    u = u + bf(c)*exp(-a*p*c^2).*sin(c*x);
    end
    plot(x,u)
    hold on
    drawnow
    xlabel('coordinata spaziale')
    ylabel('temperatura')
    title('Andamento di T nel tempo')
end










