function [ riv,outlets,shp,mat,IC ] = read_riv( riv_file_and_path )

if ( nargin~=1 )
    fprintf('[riv,outlets,shp,mat,IC] = readriv (pwd)\n\n');
    fprintf('Follow the usage above,sweet heart! :-)\n\n\n');
    return
end

fid=fopen(riv_file_and_path,'r');
if (~fid)
    fprintf('File %s does NOT exist\n',fn);
    return;
end

%river
mriv=fscanf(fid,'%g',[1,1]);
riv=fscanf(fid,'%g',[11,mriv])';

%shape
shapname=fscanf(fid,'%s',[1,1]);
mshp=fscanf(fid,'%g',[1,1]);
shp=fscanf(fid,'%g',[4,mshp])';

%Material
matname=fscanf(fid,'%s',[1,1]);
mmat=fscanf(fid,'%g',[1,1]);
mat=fscanf(fid,'%g',[6,mmat])';

%Material
ICname=fscanf(fid,'%s',[1,1]);
mIC=fscanf(fid,'%g',[1,1]);
IC=fscanf(fid,'%g',[2,mIC])';

fclose(fid);

outlets=riv(find(riv(:,4)<0),1);

end

