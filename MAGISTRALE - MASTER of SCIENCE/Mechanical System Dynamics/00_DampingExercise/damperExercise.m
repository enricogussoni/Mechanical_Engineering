clear all
%% 
%  
% 
%  
% 
% We consider a 30 metres long IPE200 beam
% 
% We define the physical properties of the beam:

L=8;          %m
E=200e9;       %N/m^2
I=1.9430e-05;   %m^4
m=22.4;        %Kg/m

EI=E*I;
%% 
% We investigate the frequency range from 0 to 20 Hz

f=linspace(0,300,100000);
omega=2*pi*f;

%% 
% $$\gamma=\sqrt[4]{\frac{m}{EI}}\sqrt{\omega}$$

gamma=(m/EI)^(1/4)*omega.^(1/2);
%% 
% We consider the usual domain equation:
% 
% $$w=\left[A \sin \left( \gamma x \right) +B\cos\left( \gamma x \right)+C\sinh\left( 
% \gamma x \right)+D\cosh\left( \gamma x \right)\right]\cos\left(\omega t + \psi 
% \right)$$
% 
% We need to impose the boundary conditions:
% 
% Null dispacement in the first pin
% 
% $$w|_0=0 \; \; \Rightarrow \;\;B+D=0$$
% 
% Null bending moment in the first pin
% 
% $$EI \frac{\partial^2w}{\partial x^2}|_0=0 \; \; \Rightarrow \;\;EI(-\gamma^2B+\gamma^2 
% D)=0$$
% 
% $$\Rightarrow \;\;-B+ D=0$$
% 
% Null dispacement in the second pin
% 
% $$w|_L=0 \;\;\Rightarrow\;\;A \sin ( \gamma L ) +B\cos( \gamma L )+C\sinh( 
% \gamma L )+D\cosh( \gamma L ) =0$$
% 
% Null bending moment in the second pin
% 
% $$EI w|_L=0 \;\;\Rightarrow\;\;-\gamma^2A \sin ( \gamma L ) -\gamma^2B\cos( 
% \gamma L )+\gamma^2C\sinh( \gamma L )+\gamma^2D\cosh( \gamma L ) =0$$
% 
% $$\Rightarrow\;\;-A \sin ( \gamma L ) -B\cos( \gamma L )+C\sinh( \gamma 
% L )+D\cosh( \gamma L ) =0$$
% 
% We can write the system BC matrix:
% 
% $$\pmatrix{0 & 1 & 0 & 1 \cr 0 & -1 & 0 & 1 \cr \sin(\gamma L) & \cos(\gamma 
% L) & \sinh(\gamma L) & \cosh(\gamma L) \cr -\sin(\gamma L) & -\cos(\gamma L) 
% & \sinh(\gamma L) & \cosh(\gamma L)}$$ 

H=@(gamma) [  0             1            0             1            ;
    0            -1            0             1            ;
    sin(gamma*L)  cos(gamma*L) sinh(gamma*L) cosh(gamma*L);
    -sin(gamma*L) -cos(gamma*L) sinh(gamma*L) cosh(gamma*L)];
%% 
% We obtain a 4 x 4 x 50 matrix containing the matrix H for each value of 
% gamma.
% 
% We are looking for the value of gama for which the determinant of H is 
% equal to zero

for i=1:length(gamma);
    dets(i)=det(H(gamma(i)));
end

semilogy(f,abs(dets),'-b')
hold on, xlabel('f [Hz]')
%% 
% I have five values of f for which the determinant is (close) to zero

i_nat=[];
for i=2:length(dets)-1
    if abs(dets(i)) < abs(dets(i-1)) && abs(dets(i)) < abs(dets(i+1))
        i_nat(end+1)=i;
    end
end

fprintf('Natural frequencies:\n %f\n%f\n%f\n%f\n%f\n',f(i_nat))
fprintf('Gamma:\n %f\n%f\n%f\n%f\n%f\n',gamma(i_nat))
plot(f(i_nat),abs(dets(i_nat)),'or')
%% 
% Now that we know the value for which the system is singular (i.e. admits 
% non-trivial solutions), we can find the modal shapes solving the *reduced system:*
% 
% $$\left[ \matrix{-&- \cr\hat{E_i} & \hat{H}(\gamma_i)&}\right]\left\{\matrix{1 
% \cr \hat{C_i} } \right\}= 0$$
% 
% $$\hat{C_i}= -\hat{H}^{-1}(\gamma_i) \hat{E_i}$$

for i_mode=1:length(i_nat)
    fprintf('MODO %i',i_mode)
    gamma_i=gamma(i_nat(i_mode));
    Hi=H(gamma_i)
    Hi_hat=Hi(2:4,2:4)
    Ei_hat=Hi(2:4,1)
    Ci_hat=[1; -Hi_hat\Ei_hat]
    
    C_hat(:,i_mode)=Ci_hat;
end
C_hat
%% 
% We can now compute the mode shapes
% 
% $$\phi_i(x)=\left[A \sin \left( \gamma_i x \right) +B\cos\left( \gamma_i 
% x \right)+C\sinh\left( \gamma_i x \right)+D\cosh\left( \gamma_i x \right)\right]$$


x=linspace(0,L,1000);
dx=x(2);
for i_mode=1:length(i_nat)
    gamma_i=gamma(i_nat(i_mode));
    phi(i_mode,:)= C_hat(1,i_mode)*sin(gamma_i*x)  + C_hat(2,i_mode)*cos(gamma_i*x) +...
        C_hat(3,i_mode)*sinh(gamma_i*x) + C_hat(4,i_mode)*cosh(gamma_i*x);
end

figure
plot(x,phi)
hold on, plot([0 L],[0 0],'--k')
ylim([-1 1])
legend({'Mode 1','Mode 2','Mode 3','Mode 4','Mode 5'},'Location','NorthOutside')
%% 
% Now we can compute the energy functions using the Lagrange equation:
% 
% $$E_K= \frac{1}{2} \int_0^L m \frac{\partial w}{\partial t}\frac{\partial 
% w}{\partial t} dx=$$
% 
% $$=\frac{1}{2} \dot{q}^T \int_0^L m \phi(x) \phi^T(x) \,dx \;\dot{q}=$$
% 
% $$=\frac{1}{2} \dot{q}^T [M] \;\dot{q}$$
% 
% $$\Rightarrow [M] = \int_0^L m \phi(x) \phi^T(x) \,dx$$

M=m*phi*phi'*dx
%% 
% As you can see *M is diagonal*
% 
% $$V= \frac{1}{2} \int_0^L EJ \frac{\partial^2 w}{\partial x^2}\frac{\partial^2 
% w}{\partial x^2} dx=$$
% 
% $$=\frac{1}{2} q^T \int_0^L EJ \phi''(x) \phi''^T(x) \,dx \;q=$$
% 
% $$=\frac{1}{2} q^T [K] \;q$$
% 
% $$\Rightarrow [K] = \int_0^L EJ \phi''(x) \phi''^T(x) \,dx$$

phi2=diff(phi,2,2)/(dx^2);

K=EI*phi2*phi2'*dx
%% 
% *As you can see K is diagonal too!*
% 
% **
% 
% *The equation of motion for this system is therefore:*
% 
% $$[M]\ddot{q} + [K]\dot{q} = 0$$
% 
% Now consider the same beam, this time with a damper in the middle:
% 
%  
% 
% We need to compute the damping function.
% 
% We say that:
% 
% $$D =\frac{1}{2} r \frac{\partial w}{\partial t}|_{x=2} \frac{\partial 
% w}{\partial t}|_{x=2}=$$
% 
% $$= \frac{1}{2} \dot{q}^T r \phi(x=2)\phi^T(x=2)\dot{q}=$$
% 
% $$= \frac{1}{2} \;\dot{q}^T [C]\;\dot{q}$$
% 
% $$\Rightarrow\; [C]=r\,\phi(x=2)\,\phi^T(x=2)$$

xDamper=2;
r=.10 * 2 * M(1,1) * omega(i_nat(1));

[~,i_damper]=min(abs(x-xDamper));
phi_damper=phi(:,i_damper)

plot(xDamper,phi_damper,'or')

R=r*phi_damper*phi_damper'*dx
%% 
% *R IS NOT DIAGONAL!!*
% 
% We want to compute the eigenvalues and eigenvectors of the damped system.
% 
% We use the the system state variable
% 
% $$z=\left\{ \matrix{ \dot{q} \cr q } \right\}$$
% 
% to write the system as
% 
% $$\matrix{[M] \ddot{q} & +  [R]\dot{q} & +  [K]q &=0 \cr[M]\dot{q} && -[M]\dot{q} 
% &=0}$$
% 
% that allow the system to be written as
% 
% $$\left[ \matrix{ [M] & [0] \cr[0] & [M]} \right]\left\{ \matrix{\ddot{q} 
% \cr\dot{q}} \right\}+\left[ \matrix{ [R] & [K] \cr-[M] & [0]} \right]\left\{ 
% \matrix{\dot{q} \crq} \right\}=\left\{ \matrix{0 \cr0} \right\}$$
% 
% $$[B]\dot{z} + [C]z=0$$
% 
% $$[A] = -[B]^{-1}[C]$$
% 
% $$\dot{z}-[A]z=0$$

B=[M,              zeros(size(M));
    zeros(size(M)), M             ];
C=[ R, K             ;
    -M, zeros(size(M))];
A=-inv(B)*C;
%% 
% we can use this equation to write the eigenvalue problem
% 
% $$[[A] -\lambda[I]]Z=0$$


[phi_damped,lamba_damped]=eig(A,'vector')
%% 
% We have 10 eigen-vectors and 10 eigen-values complex coniugated two by 
% two. We keep only one of each

phi_damped=phi_damped(:,1:2:end);
lamba_damped=lamba_damped(1:2:end);
%% 
% Since we are interested only in the first 5 lines of the previous system, 
% we keep only the first 5 values of each eigen-vector

phi_damped=phi_damped(1:5,:);
%% 
% Now we sort the eigen-values smallest-to-biggest and the eigen-vectors 
% accordingly.

[lamba_damped,sortInd]=sort(lamba_damped);
phi_damped=phi_damped(:,sortInd)
%% 
% We can plot these eigen-vetors in the complex plane:

for i=1:5
    figure
    quiver(zeros(5,1),zeros(5,1),real(phi_damped(:,i)),imag(phi_damped(:,i)))
    title(sprintf('Mode %i',i))
    xlim([-1 1]), ylim([-1 1])
    xlabel('Re(\phi)'), ylabel('Im(\phi)')
    grid on
end
%% 
% Consider Mode 2:
% 
% The damped mode will be:
% 
% $$w_2(x,t)=\phi(x)^T \phi_{D,2} \;q(t) = \phi(x)^T \phi_{D,2} \;e^{i\omega_1 
% t}$$
% 
% since $\phi_{D,2}$ is complex, also the vibrating mode will be complex, 
% with some components out of phase with respect to the "pure" mode. To look at 
% it we need to see it in the time domain:

mode=2;

figure, hold on
title('2nd mode')
ylim([-1 1])
plot(x,phi(mode,:),':k','DisplayName','Un-damped system')
plot(x,zeros(size(x)),'-k','HandleVisibility','off')
h=plot(x,zeros(size(x)),'DisplayName','Damped system');
legend show

or t=linspace(0,5/f(i_nat(mode)),500)
    if ishandle(h)
        w=real(phi'*phi_damped(:,mode)*exp(1i*omega(i_nat(mode))*t));
        h.YData=w;
        pause(.05)
    else
        return
    end
end