function [ g ] = read_geol( geol_file_and_path, geol  )

if ( nargin<1 )
    fprintf('[g]=rwgeol(cdir,[geo])\n');
    fprintf('\t[geo] is only used for writing\n\n\n');
    return
end

if (nargin >1)
    wif=1;
else
    wif=0;
end

if (~exist(geol_file_and_path,'file'))
    fprintf('File %s does NOT exist\n',fn);
    return;
end

if (wif)
    [m,n]=size(geol);
    fid=fopen(fn,'w');
    fprintf(fid,'%g\n',m);
    for i=1:m
        fprintf(fid,'%10g\t',geol(i,:));
        fprintf(fid,'\n');
    end
    g=geol;
else
    fid=fopen(geol_file_and_path,'r');
    a=fscanf(fid,'%g');
    n=a(1);
    g=reshape(a(2:end),(length(a)-1)/n,n)';
end

fclose(fid);
end

