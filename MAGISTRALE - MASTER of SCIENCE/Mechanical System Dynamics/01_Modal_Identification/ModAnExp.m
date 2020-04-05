close all;
clear all;

%% Natural frequencies
load('DH1.mat')
% figure(1)
% subplot(3,2,1); plot(freq,abs(frf(:,1)));
% subplot(3,2,3); plot(freq,abs(frf(:,2)));
% subplot(3,2,5); plot(freq,abs(frf(:,3)));
% 
% subplot(3,2,2); plot(freq,angle(frf(:,1)));
% subplot(3,2,4); plot(freq,angle(frf(:,2)));
% subplot(3,2,6); plot(freq,angle(frf(:,3)));

% load('DH2.mat')
% figure(2)
% subplot(3,1,1); plot(freq,frf(:,1));
% subplot(3,1,2); plot(freq,frf(:,2));
% subplot(3,1,3); plot(freq,frf(:,3));
%
% load('RH1.mat')
% figure(3)
% subplot(3,1,1); plot(freq,frf(:,1));
% subplot(3,1,2); plot(freq,frf(:,2));
% subplot(3,1,3); plot(freq,frf(:,3));

ch1=abs(frf(:,1));
ch2=abs(frf(:,2));
ch3=abs(frf(:,3));

[picchi1,pos1]=findpeaks(ch1,'minpeakheight',100,'minpeakdistance',300);
[picchi2,pos2]=findpeaks(ch2,'minpeakheight',100,'minpeakdistance',300);
[picchi3,pos3]=findpeaks(ch3,'minpeakheight',100,'minpeakdistance',300);

% figure(2); hold on;
% plot(freq,ch1);
% plot(freq,ch2);
% plot(freq,ch3);
% plot(freq(pos1),picchi1,'ro')
% plot(freq(pos2),picchi2,'ro')
% plot(freq(pos3),picchi3,'ro')

[picchi,pos]=findpeaks(ch1+ch2+ch3,'minpeakheight',100,'minpeakdistance',200);
nat_freq=freq(pos)

% %% DAMPING RATIO
% % Half-Power points
% sp=500;
% 
% G=sqrt(2)/2*ch1(pos);
% for ii=1:length(pos)
%     
%     act=pos(ii);
%     val1=(find((ch1(act-sp:act+sp))>G(ii),1,'first')+act-sp-1);
%     val2=(find((ch1(act-sp:act+sp))>G(ii),1,'last')+act-sp-1);
%     om1=(freq(val1)-freq(val1-1))/(ch1(val1)-ch1(val1-1))*(G(ii)-ch1(val1-1))+freq(val1-1);
%     om2=(freq(val2+1)-freq(val2))/(ch1(val2+1)-ch1(val2))*(G(ii)-ch1(val2))+freq(val2);
%     psiHP(ii,1)=(om2^2-om1^2)/(4*nat_freq(ii)^2);
% end
% G=sqrt(2)/2*ch2(pos);
% for ii=1:length(pos)
%     
%     act=pos(ii);
%     val1=(find((ch2(act-sp:act+sp))>G(ii),1,'first')+act-sp-1);
%     val2=(find((ch2(act-sp:act+sp))>G(ii),1,'last')+act-sp-1);
%     om1=(freq(val1)-freq(val1-1))/(ch2(val1)-ch2(val1-1))*(G(ii)-ch2(val1-1))+freq(val1-1);
%     om2=(freq(val2+1)-freq(val2))/(ch2(val2+1)-ch2(val2))*(G(ii)-ch2(val2))+freq(val2);
%     psiHP(ii,2)=(om2^2-om1^2)/(4*nat_freq(ii)^2);
% end
% G=sqrt(2)/2*ch3(pos);
% for ii=1:length(pos)
%     
%     act=pos(ii);
%     val1=(find((ch3(act-sp:act+sp))>G(ii),1,'first')+act-sp-1);
%     val2=(find((ch3(act-sp:act+sp))>G(ii),1,'last')+act-sp-1);
%     om1=(freq(val1)-freq(val1-1))/(ch3(val1)-ch3(val1-1))*(G(ii)-ch3(val1-1))+freq(val1-1);
%     om2=(freq(val2+1)-freq(val2))/(ch3(val2+1)-ch3(val2))*(G(ii)-ch3(val2))+freq(val2);
%     psiHP(ii,3)=(om2^2-om1^2)/(4*nat_freq(ii)^2);
% end
% format shortEng
% psiHP
% 
% % Slope of the phase
% slope=(angle(frf(pos+1,1))-angle(frf(pos-1,1)))/(0.04);
% psiSP(:,1)=-1./(nat_freq.*slope);
% slope=(angle(frf(pos+1,2))-angle(frf(pos-1,2)))/(0.04);
% psiSP(:,2)=-1./(nat_freq.*slope);
% slope=(angle(frf(pos+1,3))-angle(frf(pos-1,3)))/(0.04);
% psiSP(:,3)=-1./(nat_freq.*slope);
% psiSP


%% MODE SHAPES
figure(3)
hold on
load('modes.mat')
% plot(modes.x,modes.y)
x1=0.105;
x2=0.415;
x3=0.6;
L=1.2;
xdis=linspace(0,L,1000);
undeformed=@(x) x*0;

pti(:,1)=ch1(pos)./angle(frf(pos,1));
pti(:,2)=ch2(pos)./angle(frf(pos,2));
pti(:,3)=ch3(pos)./angle(frf(pos,3));
%pti1=pti1./(max(abs(pti1)*max(modes.y)));
%pti2=pti2./(max(abs(pti2)));
%pti3=pti3./(max(abs(pti3)));

for ii=1:length(nat_freq)
    switch ii
        case 1
            pti(ii,:)=-pti(ii,:)/max(abs(pti(ii,:)))*max(modes.y(:,ii));
        case 2
            pti(ii,:)=pti(ii,:)/max(abs(pti(ii,:)))*max(modes.y(:,ii));
        case 3
            pti(ii,:)=-pti(ii,:)/max(abs(pti(ii,:)))*max(modes.y(:,ii));
        case 4
            pti(ii,:)=-pti(ii,:)/max(abs(pti(ii,:)))*max(modes.y(:,ii));
        case 5
            pti(ii,:)=pti(ii,:)/max(abs(pti(ii,:)))*max(modes.y(:,ii));
    end
end

for ii=1:5
    figure(2+ii);
    hold on;
    plot(xdis,undeformed(xdis),'--k','linewidth',2);
    plot(modes.x,modes.y(:,ii));
    plot(x1,pti(ii,1),'ro')
    plot(x2,pti(ii,2),'go')
    plot(x3,pti(ii,3),'yo')
    legend('Analytical','Experimental');
    
end
%% CHECK
figure(10)
plot(freq,cohe(:,1))
figure(11)
plot(freq,cohe(:,2))
figure(12)
plot(freq,cohe(:,3))


