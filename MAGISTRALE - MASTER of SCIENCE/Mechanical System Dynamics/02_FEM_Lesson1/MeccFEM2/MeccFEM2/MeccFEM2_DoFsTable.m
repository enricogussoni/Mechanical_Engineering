function MeccFEM2_DoFsTable(idb)
%MECCFEM2_DOFSTABLE Display the dofs table of the system
%   MECCFEM2_DOFSTABLE(IDB) displais the DoFs table reading the 3-by-N idb 
%   matrix containing the global dof index for each node's dofs.


figure
data=[(1:size(idb,1))' idb];
cnames = {'Node','X DoF','Y DoF','Rotation'};
uitable('Data',data,...
          'ColumnName',cnames,...
          'RowName',{},...
          'Units','normalized',...
          'Position',[.05 .05 .9 .9]);
       