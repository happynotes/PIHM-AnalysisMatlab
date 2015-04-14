function [  msh,pt ] = read_mesh( mesh_file_and_path )

if ( nargin~=1 )
    fprintf('Usage:\n \t[ mesh,pt ] = readmesh (mesh_file_and_path)\n\n');
    return
end


fid=fopen(mesh_file_and_path,'r');

if (fid<0)
    fprintf('File %s does NOT exist\n',fn);
    return;
end

a=fscanf(fid,'%g');

if( isempty(a))
    fprintf('Error:File .mesh is empty\n\n');
    return
end

mm=a(1);
mp=a(2);

%Why use these sizes? Specific reason?
nm=7;
np=5;

%http://www.mathworks.com/help/matlab/ref/reshape.html?searchHighlight=reshape
msh=reshape(a(3:mm*nm+2),nm,mm)';
pt=reshape(a(mm*nm+3:end),np,mp)';

fclose(fid);

end

