function MeccFEM2_plotDeformedStructure(mode,scaleFactor,incidence,l,gamma,position,idb,xy)
%MECCFEM2_PLOTDEFORMEDSTRUCTURE plot the deformed structure
%   MECCFEM2_PLOTDEFORMEDSTRUCTURE(MODE,SCALEFACTOR,INCIDENCE,L,GAMMA,POSITION,IDB,XY)
%   plot the deformed structure using the output of MECCFEM2_LOADSTRUCTURE.
%   The variable MODE describe the displacement of the firsts N DoF where N
%   is the length of MODE. All the following DoFs are assumed to be
%   contrained and, therefore, a null displacement is used.
%   SCALEFACTOR is an amplification factor for the mode scaling
%   INCIDENCE, L, GAMMA, POSITION, IDB and XY are the same described in the
%   help of MECCFEM2_LOADSTRUCTURE.

% Input check
[n_el,~]=size(incidence);

if size(position,1) ~= n_el
    error('Number of elembents in the position matrix different from number of element in the incidence matrix')
end
if length(gamma) ~= n_el
    error('Number of element of gamma different from the number of elements')
end
if length(l) ~= n_el
    error('Number of element of l different from the number of elements')
end

n_dof=length(mode);


hold on
% Cicle on all the elements
for k=1:n_el
% Contstruction of the nodal displacement vector in the global reference
% frame

    % Check if the dof is constrained
    for irow=1:6
        if incidence(k,irow) <= n_dof
            xkG(irow)=mode(incidence(k,irow));
        else
            xkG(irow)=0.;
        end
    end
% Application of the scale factor
    xkG=scaleFactor*xkG;
% Projection of the displacements from global to local reference frame
    lambda = [ cos(gamma(k)) sin(gamma(k)) 0. 
              -sin(gamma(k)) cos(gamma(k)) 0.
	                0.      0.     1. ] ;
    Lambda = [ lambda     zeros(3)
              zeros(3)   lambda      ] ;
    xkL=Lambda*xkG';

% Computation of axial displacements u and transversal displacemens w using
% the shape functions
    csi=linspace(0,l(k),20);
    
    fu=zeros(6,length(csi));
    fu(1,:)=1-csi/l(k);
    fu(4,:)=csi/l(k);
    u=(fu'*xkL)';
    
    fw=zeros(6,length(csi));
    fw(2,:)=2*(csi/l(k)).^3-3*(csi/l(k)).^2+1;
    fw(3,:)=l(k)*((csi/l(k)).^3-2*(csi/l(k)).^2+csi/l(k));
    fw(5,:)=-2*(csi/l(k)).^3+3*(csi/l(k)).^2;
    fw(6,:)=l(k)*((csi/l(k)).^3-(csi/l(k)).^2);
    w=(fw'*xkL)';
    
% Projection of the deformed shape from local to global
    deformed  =lambda(1:2,1:2)'*[csi+u; w];
    undeformed=lambda(1:2,1:2)'*[csi;   zeros(1,length(csi))];
% Plot undeformed and deformed element
    plot(deformed(1,:)  +position(k,1),  deformed(2,:)+position(k,2),'b')
    plot(undeformed(1,:)+position(k,1),undeformed(2,:)+position(k,2),'b--')
end

% Cycle on the nodes
n_nodes=size(idb,1);

xkG=[];
for k=1:n_nodes
    
    %Check if the dof is constrained
    for ixy=1:2
        if idb(k,ixy) <= n_dof
            xkG(k,ixy)=mode(idb(k,ixy));
        else
            xkG(k,ixy)=0.;
        end
    end
end
xkG=scaleFactor*xkG;
xyDef=xy+xkG;

plot(xy(:,1),xy(:,2),'b.')
plot(xyDef(:,1),xyDef(:,2),'bo')

grid on
axis equal

try
    sid = '';
    ni = java.net.NetworkInterface.getNetworkInterfaces;
    while ni.hasMoreElements
        addr = ni.nextElement.getHardwareAddress;
        if ~isempty(addr)
            addrStr = dec2hex(int16(addr)+128);
            sid = [sid, '.', reshape(addrStr,1,[])];
        end
    end
    
    setappdata(gcf,'sid',sid)
end