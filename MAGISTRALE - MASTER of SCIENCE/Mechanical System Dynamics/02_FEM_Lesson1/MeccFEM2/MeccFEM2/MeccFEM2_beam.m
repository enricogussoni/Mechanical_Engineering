function [mG,kG] = MeccFEM2_beam (l,m,EA,EJ,gamma)
%MECCFEM2_BEAM Bernulli's beam stiffness and mass matrix
%   [MG,KG] = MECCFEM2_BEAM (L,M,EA,EJ,GAMMA) generates the stiffness
%   matrix KG and the mass matrix MG in the global reference frame for a
%   bernulli beam element with length L, mass per unit length M, axial
%   stiffness EA, bending stiffness EJ and rotated by an angle GAMMA
%   (defined in radiants counterclock-wise)

%% Mass matrix
% Mass matrix in local reference frame
mL = m*l*[1./3.  0.          0.          1./6.  0.          0.
          0.     13./35.     11.*l/210.  0.     9./70.      -13*l/420.
          0.     11.*l/210.  l^2/105.    0.     13*l/420.   -l^2/140.
          1./6.  0.          0.          1./3.  0.          0.
          0.     9./70.      13*l/420.   0.     13./35.     -11.*l/210.
          0.     -13*l/420.  -l^2/140.   0.     -11.*l/210. l^2/105.   ] ;

%% Stiffness matrix
% Stiffness matrix in local reference frame
% contribution of the axial stiffness
kL_ax = EA/l* [ 1 0 0 -1 0 0
                0 0 0  0 0 0 
                0 0 0  0 0 0 
                -1 0 0  1 0 0 
                0 0 0  0 0 0 
                0 0 0  0 0 0 ] ; 

% contribution of the bending stiffness
kL_fl = EJ * [ 0.    0.       0.      0.    0.       0.     
               0.  +12./l^3  6./l^2   0.  -12./l^3  6./l^2
               0.   6./l^2  +4./l     0.   -6./l^2  +2./l
               0.    0.       0.      0.    0.       0. 
               0.  -12./l^3  -6./l^2  0.  +12./l^3  -6./l^2
               0.   6./l^2  +2./l    0.   -6./l^2  +4./l    ] ;
           
kL = kL_ax+kL_fl ;

%% Transformation into the global reference frame

% 3x3 lambda matrix
lambda = [ cos(gamma) sin(gamma) 0. 
          -sin(gamma) cos(gamma) 0.
              0.         0.      1.] ;

% 6x6 lambda matrix for the complete element
Lambda = [ lambda     zeros(3)
           zeros(3)   lambda      ] ;
       
 % the local reference sistem is rotated while the global remains fixed
mG = Lambda' * mL * Lambda ; 
kG = Lambda' * kL * Lambda ;