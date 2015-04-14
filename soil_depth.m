function [ sd,eid ] = soil_depth( mesh_file_and_path )

if ( nargin < 1 )
    fprintf('Usage:\n \t [sd,eid]=sdepth(mesh_file_and_path)\n\n');
    return
end


[msh,pt]=read_mesh(mesh_file_and_path);
m=size(msh,1);
%elv=(pt(msh(:,2),5)+pt(msh(:,3),5)+ pt(msh(:,4),5) )/3;
elv=zeros(m,1);
belv=zeros(m,1);
for i=1:3
    elv=elv+pt(msh(:,i+1),5)./3;
end
for i=1:3
    belv=belv+pt(msh(:,i+1),4)./3;
end

sd=elv-belv;

eid=find(sd<=0);
if ~isempty(eid)
    disp('Cells have zero/negative soil depth');
    fprintf('\n%d cells have zero/negative soil depth\n\n', eid);
end

end

