function [file_i,xy,nnod,sizew,idf,ndof,incidence,l,gamma,m,EA,EJ,position,nbeam]=MeccFEM2_loadstructure(file_i)
%MECCFEM2_LOADSTRUCTURE load the structure data
%   [FILE_NAME,XY,NNOD,SIZEW,IDB,NDOF,INCIDENCE,L,GAMMA,M,EA,EJ,POSIZ,NBEAM]=MECCFEM2_LOADSTRUCTURE
%   loads the structure information from a given file. The outputs are:
%   
%   * FILE_NAME is the input file name
%   * XY are the coordinates of the nodes
%   * NNOD is the number of nodes
%   * SIZEW is the diametre of the circumscribed circle to the structure
%     (used for plotting)
%   * IDF is a matrix containing the index of the 3 degrees of freedom of
%     each node
%   * NDOF is the number of dofs in the structure
%   * INCIDENCE is the incidence matrix for the system
%   * L GAMMA, M, EA, EJ are vectors containing the length, the rotation, the 
%     mass per unit length, the axial stiffness and the bending stiffness of 
%     each element
%   * POSITION is a 4 column matrix containing on the n-th row the x,y
%     coordinates of the first node of the n-th element and the x,y
%     coordinated of the last node of the n-th element
%   * NBEAM is the number of beams of the structure
%
%   [...]=MECCFEM2_LOADSTRUCTURE(FILENAME) loads the data from file
%   FILENAME without prompting for the file name


%#ok<*AGROW>

%% File loading
if nargin < 1
    file_i=input('\nStructure input file''s name (without .inp) = ','s') ;
else
    if length(file_i)>4 && strcmp(file_i(end-4:end),'.inp')
        file_i=file_i(1:end-5);
    end
end

% Check for file existance
if ~exist(sprintf('%s.inp',file_i),'file')
    error('File %s.inp does not exists',file_i)
end

% Open the file
fid=fopen(sprintf('%s.inp',file_i));

%% NODES reading

%Skip at the beginning of the NODES section
MeccFEM2_findcard(fid,'*NODES')

nnod=0;
while true    
    
    %Get the first un-commented line
    line=MeccFEM2_fgetul(fid);
    
    if isempty(strfind(line,'*ENDNODES'))
        %Extract the line's data
        line_data=sscanf(line,'%i %i %i %i %f %f')';
        % Increase the counter
        nnod=nnod+1;
        
        %check for consistency
        if length(line_data)<6
            error('Invalid data format in NODES section, node number %i, column %i',nnod,length(line_data))
        end
        if nnod ~=line_data(1)
            error('Nodes are not numbered progressively')
        end
        
        %Extract contrain information
        iconst(nnod,:)=line_data(2:4); 
        
        %Extract coordinates
        xy(nnod,:)=line_data(5:6);
    else
        break
    end
end

fprintf('Number of nodes in the structure:\t%i\n',nnod)

sizew=sqrt((max(xy(:,1))-min(xy(:,1)))^2+(max(xy(:,2))-min(xy(:,2)))^2);

%% CONTRUCTION OF IDB MATRIX AND DOF NUMERATION

% Construction of IDB matrix of the nodes numeration
% 1) Free dof

ndof=0;

%cycling on each node
for i=1:nnod
    %cycling on each coordinate
    for j=1:3
        %in the dof is free (0) call it with the following name
        if iconst(i,j) == 0
            ndof=ndof+1;
            idf(i,j)=ndof;
        end
    end
end

fprintf('Number of d.o.f of the structure:\t%i\n',ndof)

% 2) constrained dof

ndoc=ndof;
%cycling on each node
for i=1:nnod
    %cycling on each coordinate
    for j=1:3
        %in the dof is constrained (1) call it with the following name
        if iconst(i,j) == 1
            ndoc=ndoc+1;
            idf(i,j)=ndoc;
        end
    end
end


%% BEAMS READING
MeccFEM2_findcard(fid,'*BEAMS')
nbeam=0;
while true
    %Get the first un-commented line
    line=MeccFEM2_fgetul(fid);
    
    if isempty(strfind(line,'*ENDBEAMS'))
        %Extract the line's data
        line_data=sscanf(line,'%i %i %i %f %f %f')';
        
        % Increase the counter
        nbeam=nbeam+1;
        
        %check for consistency
        if length(line_data)<6
            error('Invalid data format in BEAMS section, beam number %i, column %i',nbeam,length(line_data))
        end
        if nbeam ~=line_data(1)
            error('Beams are not numbered progressively')
        end
        
        %Extract the initial and final nodes
        nodes(nbeam,:)=line_data(2:3);
        
        %Extract fisical variables
        m(1,nbeam)=line_data(4);
        EA(1,nbeam)=line_data(5);
        EJ(1,nbeam)=line_data(6);
        
        %Extract initial and final coordinates
        position(nbeam,:)=[xy(nodes(nbeam,1),:) xy(nodes(nbeam,2),:)];
        deltaX=diff(position(nbeam,[1 3]));
        deltaY=diff(position(nbeam,[2 4]));
        
        %Compute the length with pitagora
        l(1,nbeam)=sqrt(deltaX^2+deltaY^2);
        
        %Compute the angle with atan2
        gamma(1,nbeam)=atan2(deltaY,deltaX);
        
        %Fill the incidence matrix using the dof indexes of the first and
        %the second node
        incidence(nbeam,:)=[idf(nodes(nbeam,1),:) idf(nodes(nbeam,2),:)];
    else
        break
    end
end

fprintf('Number of beams in the structure:\t%i\n',nbeam)

%% CLOSING
% Close the input file.
fclose(fid) ;