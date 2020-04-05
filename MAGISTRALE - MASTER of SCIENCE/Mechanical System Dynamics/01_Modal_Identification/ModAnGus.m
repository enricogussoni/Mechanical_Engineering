%% Fundamentals of experimental modal analysis
% Mechanical Systems Dynamics

% DH1.mat hammer in DH1 position
% DH2.mat hammer in DH2 position
% RDH1.mat Acc1 and DH1 position interchanged
% (reciprocity against DH1.mat = If I swap the position of the imput force 
% and the measurement point the FRF remains the same.)
%
% freq: frequency vector (resolution 0.02 Hz)
% frf: frequency response functions (complex), collected by
% columns (A1,A2,A3)
% cohe: coherence functions, collected by columns (A1,A2,A3)

% Accelerometers are in A1, A2, A3 and hammering in DH1 and DH2

% ______A1______A2______A3_____DH1______DH2
%        |   ^   |       |      |   ^    |
% <-----> 105 mm |       |      |        |
% <--------------> 415 mm|      |        |
% <---------------------> 600 mm|        |
% <-----------------------------> 815 mm |
% <--------------------------------------> 1065 mm

% frequency resoltion is 0.02 Hz

close all;
clear all;
clc;

%% Natural frequencies
load('DH1.mat')
figure(1)
subplot(3,2,1); plot(freq,abs(frf(:,1))); xlabel('f [Hz]'); ylabel('|fft(A1)|'); % title('Accelerometer 1 - DH1');
subplot(3,2,3); plot(freq,abs(frf(:,2))); xlabel('f [Hz]'); ylabel('|fft(A2)|'); % title('Accelerometer 2 - DH1');
subplot(3,2,5); plot(freq,abs(frf(:,3))); xlabel('f [Hz]'); ylabel('|fft(A3)|'); % title('Accelerometer 3 - DH1');

subplot(3,2,2); plot(freq,angle(frf(:,1))); xlabel('< [rad]'); ylabel('<fft(A1)'); % title('Accelerometer 1 - DH1');
subplot(3,2,4); plot(freq,angle(frf(:,2))); xlabel('< [rad]'); ylabel('<fft(A2)'); % title('Accelerometer 2 - DH1');
subplot(3,2,6); plot(freq,angle(frf(:,3))); xlabel('< [rad]'); ylabel('<fft(A3)'); % title('Accelerometer 3 - DH1');

% load('DH2.mat')
% figure(2)
% subplot(3,1,1); plot(freq,frf(:,1)); xlabel('f [Hz]'); ylabel('|fft(A1)|'); % title('Accelerometer 1 - DH2'); 
% subplot(3,1,2); plot(freq,frf(:,2)); xlabel('f [Hz]'); ylabel('|fft(A2)|'); % title('Accelerometer 2 - DH2');
% subplot(3,1,3); plot(freq,frf(:,3)); xlabel('f [Hz]'); ylabel('|fft(A3)|'); % title('Accelerometer 3 - DH2');
%
% subplot(3,2,2); plot(freq,angle(frf(:,1)));  xlabel('< [rad]'); ylabel('<fft(A1)'); % title('Accelerometer 1 - DH2');
% subplot(3,2,4); plot(freq,angle(frf(:,2)));  xlabel('< [rad]'); ylabel('<fft(A2)'); % title('Accelerometer 2 - DH2');
% subplot(3,2,6); plot(freq,angle(frf(:,3)));  xlabel('< [rad]'); ylabel('<fft(A3)'); % title('Accelerometer 3 - DH2');
%
% load('RH1.mat')
% figure(3)
% subplot(3,1,1); plot(freq,frf(:,1)); xlabel('f [Hz]'); ylabel('|fft(A1)|'); % title('Accelerometer 1 - RH1');
% subplot(3,1,2); plot(freq,frf(:,2)); xlabel('f [Hz]'); ylabel('|fft(A1)|'); % title('Accelerometer 1 - RH1');
% subplot(3,1,3); plot(freq,frf(:,3)); xlabel('f [Hz]'); ylabel('|fft(A1)|'); % title('Accelerometer 1 - RH1');
%
% subplot(3,2,2); plot(freq,angle(frf(:,1))); xlabel('< [rad]'); ylabel('<fft(A1)'); % title('Accelerometer 1 - RH1');
% subplot(3,2,4); plot(freq,angle(frf(:,2))); xlabel('< [rad]'); ylabel('<fft(A2)'); % title('Accelerometer 2 - RH1');
% subplot(3,2,6); plot(freq,angle(frf(:,3))); xlabel('< [rad]'); ylabel('<fft(A3)'); % title('Accelerometer 3 - RH1');

ch1=abs(frf(:,1));
ch2=abs(frf(:,2));
ch3=abs(frf(:,3));

[picchi1,pos1]=findpeaks(ch1,'minpeakheight',100,'minpeakdistance',300);
[picchi2,pos2]=findpeaks(ch2,'minpeakheight',100,'minpeakdistance',300);
[picchi3,pos3]=findpeaks(ch3,'minpeakheight',100,'minpeakdistance',300);
% [...] = findpeaks(...,'MinPeakHeight',MPH) finds only those peaks that
%     are greater than the minimum peak height, MPH. MPH is a real-valued
%     scalar. The default value of MPH is -Inf.
% [...] = findpeaks(...,'MinPeakDistance',MPD) finds peaks separated by
%     more than the minimum peak distance, MPD. This parameter may be
%     specified to ignore smaller peaks that may occur in close proximity to
%     a large local peak.

figure(2); hold on;
plot(freq,ch1);
plot(freq,ch2);
plot(freq,ch3);
plot(freq(pos1),picchi1,'ro')
plot(freq(pos2),picchi2,'ro')
plot(freq(pos3),picchi3,'ro')

[picchi,pos]=findpeaks(ch1+ch2+ch3,'minpeakheight',200,'minpeakdistance',300);
nat_freq=freq(pos)

%% DAMPING RATIO (psi)
% Half-Power points (8.2.5.3 on Advanced Dynamics of Mechanical Systems)

sp=500;               % arbitrary wide interval to define a region around 
                      % the peak to find values of ch that are > of half power level G
                      
G=sqrt(2)/2*ch1(pos); % ch1 has olready been calculated with abs(...)
for ii=1:length(pos)  % for-cycle on natural frequencies
    act=pos(ii);
    val1=(find((ch1(act-sp:act+sp))>G(ii),1,'first')+act-sp-1); % trovo due valori (1 a destra e 2 a sinistra del picco) tali che lo spetro supera il valore di G
    % 1,'first' fa si che si trovi il primo valore per cui ch1(...)>G(ii)
    val2=(find((ch1(act-sp:act+sp))>G(ii),1,'last')+act-sp-1); % la posizione viene data rispetto al picco, +act-sp-1 serve a fornire il valore effettivo (sulla scala del grafico) della frequenza 
    % 1,'last' fa si che si trovi l'ultimo valore per cui ch1(...)>G(ii)
    om1=(freq(val1)-freq(val1-1))/(ch1(val1)-ch1(val1-1))*(G(ii)-ch1(val1-1))+freq(val1-1); %???
    om2=(freq(val2+1)-freq(val2))/(ch1(val2+1)-ch1(val2))*(G(ii)-ch1(val2))+freq(val2);
    psiHP(ii,1)=(om2^2-om1^2)/(4*nat_freq(ii)^2);
end

G=sqrt(2)/2*ch2(pos);
for ii=1:length(pos)
    
    act=pos(ii);
    val1=(find((ch2(act-sp:act+sp))>G(ii),1,'first')+act-sp-1);
    val2=(find((ch2(act-sp:act+sp))>G(ii),1,'last')+act-sp-1);
    om1=(freq(val1)-freq(val1-1))/(ch2(val1)-ch2(val1-1))*(G(ii)-ch2(val1-1))+freq(val1-1);
    om2=(freq(val2+1)-freq(val2))/(ch2(val2+1)-ch2(val2))*(G(ii)-ch2(val2))+freq(val2);
    psiHP(ii,2)=(om2^2-om1^2)/(4*nat_freq(ii)^2);
end

G=sqrt(2)/2*ch3(pos);
for ii=1:length(pos)
    
    act=pos(ii);
    val1=(find((ch3(act-sp:act+sp))>G(ii),1,'first')+act-sp-1);
    val2=(find((ch3(act-sp:act+sp))>G(ii),1,'last')+act-sp-1);
    om1=(freq(val1)-freq(val1-1))/(ch3(val1)-ch3(val1-1))*(G(ii)-ch3(val1-1))+freq(val1-1);
    om2=(freq(val2+1)-freq(val2))/(ch3(val2+1)-ch3(val2))*(G(ii)-ch3(val2))+freq(val2);
    psiHP(ii,3)=(om2^2-om1^2)/(4*nat_freq(ii)^2);
end
format shortEng % Engineering format that has at least 5 digits and a power that is a multiple of three
psiHP

% Slope of the phase
slope=(angle(frf(pos+1,1))-angle(frf(pos-1,1)))/(0.04); % numerical derivate (incremental product)
psiSP(:,1)=-1./(nat_freq.*slope);
slope=(angle(frf(pos+1,2))-angle(frf(pos-1,2)))/(0.04);
psiSP(:,2)=-1./(nat_freq.*slope);
slope=(angle(frf(pos+1,3))-angle(frf(pos-1,3)))/(0.04);
psiSP(:,3)=-1./(nat_freq.*slope);
psiSP

%% MODE SHAPES
figure(3)
hold on
load('modes.mat')   % Analytical solution normalized to mode1
% plot(modes.x,modes.y)
x1=0.105;
x2=0.415;
x3=0.6;
L=1.2;
xdis=linspace(0,L,1000);
undeformed=@(x) x*0;

% For each frequency (mode shape), the ratio among the peaks |Gk(om_i)| along
% with the relative phase \(Gk(om_i) of different accelerometers k give the
% shape of the mode.
pti(:,1)=ch1(pos)./angle(frf(pos,1));
pti(:,2)=ch2(pos)./angle(frf(pos,2));
pti(:,3)=ch3(pos)./angle(frf(pos,3));
%pti1=pti1./(max(abs(pti1)*max(modes.y)));
%pti2=pti2./(max(abs(pti2)));
%pti3=pti3./(max(abs(pti3)));

% Normalization
for ii=1:length(nat_freq) 
    switch ii
        case 1 % every mode of vibration is divided for the maximum between 
               % all the modes 'max(abs(pti(ii,:))' and multiplied for it's own maximum
            pti(ii,:)=-pti(ii,:)/max(abs(pti(ii,:)))*max(modes.y(:,ii));
        case 2
            pti(ii,:)=pti(ii,:)/max(abs(pti(ii,:)))*max(modes.y(10:end-10,ii)); % 10:end-10 is needed to correct the plot after normalization
        case 3
            pti(ii,:)=-pti(ii,:)/max(abs(pti(ii,:)))*max(modes.y(10:end-10,ii));
        case 4
            pti(ii,:)=-pti(ii,:)/max(abs(pti(ii,:)))*max(modes.y(10:end-10,ii));
        case 5
            pti(ii,:)=pti(ii,:)/max(abs(pti(ii,:)))*max(modes.y(10:end-10,ii));
            
    end
end

for ii=1:length(nat_freq)
    figure(2+ii);
    hold on;
    plot(xdis,undeformed(xdis),'--k','linewidth',2);
    plot(modes.x,modes.y(:,ii));
    plot(x1,pti(ii,1),'ro')
    plot(x2,pti(ii,2),'go')
    plot(x3,pti(ii,3),'yo')
    legend('Undeformed','Analytical','Experimental');
    
end

%% CHECK

% 'cohe' tells how linear is the sistem at a specific frequency 
% (ideally is always =1). 
% cohe(f_i)<1 => the system doesn't behave linearly at that frequency.

figure(10)
hold on;
plot(freq,cohe)
legend('ch1','ch2','ch3')
