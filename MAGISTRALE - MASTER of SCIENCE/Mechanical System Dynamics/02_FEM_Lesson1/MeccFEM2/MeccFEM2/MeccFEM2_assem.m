function [M,K] = MeccFEM2_assem(incidence,l,m,EA,EJ,gamma,ndof)
%MECCFEM2_ASSEM Assemble the global FEM matrixes
%   [M,K] = MECCFEM2_ASSEM(INCIDENCE,L,M,EA,EJ,GAMMA,NDOF) will return the mass
%   and the stiffness matrixes M and K for the given system.
%   * INCIDENCE is the incidence matrix for the system
%   * L, M, EA, EJ, GAMMA are vectors containing the length, the rotation, the 
%     mass per unit length, the axial stiffness and the bending stiffness of 
%     each element
%   * NDOF is the total (free and constrained) number of dofs in the structure

[n_el,nc]=size(incidence);

if nc ~= 6
    error('Number of column in the incidence matrix different from 6')
end
if length(l) ~= n_el
    error('Number of element of l different from the number of elements')
end
if length(m) ~= n_el
    error('Number of element of m different from the number of elements')
end
if length(EA) ~= n_el
    error('Number of element of EA different from the number of elements')
end
if length(EJ) ~= n_el
    error('Number of element of EJ different from the number of elements')
end
if length(gamma) ~= n_el
    error('Number of element of gamma different from the number of elements')
end

% Warning for large elemnts number
if n_el > 50
    warning('Number of elements bigger than 50. Performances could be low')
    answer=input('Proceed anyway? (Y/N)','s');
    if strcmpi(answer,'n')
        error('Execution stopped by the user')
    end
end


% M and K assembly

M=zeros(ndof,ndof);
K=zeros(ndof,ndof);
for k=1:n_el
    % Computation of the element's K and M matrixes
    [mG,kG] = MeccFEM2_beam(l(k),m(k),EA(k),EJ(k),gamma(k));
    
    % Assembly into the global matrixes
    for irow=1:6
        for icol=1:6
            i1=incidence(k,irow); % the incidence matric tells us the correspondence
                                  % between local and global degreees of
                                  % freedom.
            i2=incidence(k,icol); % i1,i2 are i_row_glob and i_col_glob while 
                                  % irow and icol are the local indexes of
                                  % the element.
            M(i1,i2)=M(i1,i2)+mG(irow,icol);
            K(i1,i2)=K(i1,i2)+kG(irow,icol);
        end
    end
end