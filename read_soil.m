function [ ss ] = read_soil( soil_file_and_path, soil )

if ( nargin<1 )
    fprintf('[ss]=rwsoil(cdir,soil)\n\n\n');
    fprintf('Follow the usage above,sweet heart! :-)\n\n\n');
    return
end

if (nargin >1)
    wif=1;
else
    wif=0;
end

if (~exist(soil_file_and_path,'file'))
    fprintf('File %s does NOT exist\n',fn);
    return;
end

if (wif)
    [m,n]=size(soil);
    fid=fopen(soil_file_and_path,'w');
    fprintf(fid,'%g\n',m);
    for i=1:m
        fprintf(fid,'%10g\t',soil(i,:));
        fprintf(fid,'\n');
    end
    ss=soil;
else
    fid=fopen(soil_file_and_path,'r');
    a=fscanf(fid,'%g');
    n=a(1);
    ss=reshape(a(2:end),(length(a)-1)/n,n)';
end
fclose(fid);

end

