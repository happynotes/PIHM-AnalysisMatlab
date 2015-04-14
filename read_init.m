function [ minit,rinit ] = read_init( init_file_and_path ,mesh_file_and_path, minit, rinit)

if ( nargin < 2 )
    fprintf('[minit,rinit]=rwinit(init_file_and_path,mesh_file_and_path,[minit,rinit])\n\n\n');
    return
end

if (nargin > 2)
    wif=1;
else
    wif=0;
end

if (~exist(init_file_and_path,'file'))
    fprintf('File %s does NOT exist\n',fn);
    return;
end
if (~exist(mesh_file_and_path,'file'))
    fprintf('File %s does NOT exist\n',fn);
    return;
end

if (wif)
    fid=fopen(init_file_and_path,'w');
    fprintf(fid,'%10g\t%10g\t%10g\t%10g\t%10g\n',minit');
    fprintf(fid,'%10g\t%10g\n',rinit');
    fprintf('Tips: initial condition was changed.\n\t Please make sure the mesh file is match this change\n\n\n'); 
else
    msh=read_mesh(mesh_file_and_path);
    mm=size(msh,1);
    fid=fopen(init_file_and_path,'r');
    a=fscanf(fid,'%g');
    mn=5;
    minit=reshape(a(1:mm*mn),mn,mm)';
    rn=2;
    rm=(length(a)-mm*mn)/rn;
    rinit=reshape(a(mm*5+1:end),rn,rm)';
end
fclose(fid);

end

