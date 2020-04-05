function dis_stru(posiz,l,gamma,xy);
% disegna la struttura indeformata

% passo 1: disegno posizioni nodali
plot(xy(:,1),xy(:,2),'ro')
hold on

% passo 2: disegno elementi
for i=1:length(posiz)
    xin=posiz(i,1);
    yin=posiz(i,2);
    xfi=posiz(i,1)+l(i)*cos(gamma(i));
    yfi=posiz(i,2)+l(i)*sin(gamma(i));
    plot([xin xfi],[yin yfi],'b')
end
grid
title('Undeformed Structure')