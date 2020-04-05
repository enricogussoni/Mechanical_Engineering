mode=2;

figure, hold on
title('2nd mode')
ylim([-1 1])
plot(x,phi(mode,:),':k','DisplayName','Un-damped system')
plot(x,zeros(size(x)),'-k','HandleVisibility','off')
h=plot(x,zeros(size(x)),'DisplayName','Damped system');
legend show


for t=linspace(0,5/f(i_nat(mode)),500)
    if ishandle(h)
        w=real(phi'*phi_damped(:,mode)*exp(1i*omega(i_nat(mode))*t));
        h.YData=w;
        pause(.05)
    else
        return
    end
end