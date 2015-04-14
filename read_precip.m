function [ dp,dt,Pm ] = read_precip( precip_file_and_path )

maxp=100000;
if ( nargin<1 )
    fprintf('[[dp,dt,Pm]=readP(precip_file_and_path)\n\n\n');
    disp('daily precipitation(avg)');
    fprintf('\tPm\t- daily P in Maxtix');
    fprintf('Follow the usage above,sweet heart! :-)\n\n\n');
    return
end

if (~exist(precip_file_and_path,'file'))
    fprintf('File %s does NOT exist\n',fn);
    return;
end

fid=fopen(precip_file_and_path,'r');

a=fscanf(fid,'%g \n');
m=a(1);
mat=zeros(m,maxp);
for i=1:m
    head=fgets(fid);
    b=fscanf(fid,'%g %g %g %g',[4,inf]); 
    n=size(b,2);
    mat(i,1:n)=b(2,:);
end
fclose(fid);
P=mat(1:m,1:n);

p=mean(P,1);
t=b(1,1:n);

s=1;
i=2;
nt=length(t);
dp=zeros(nt,1);
dt=zeros(nt,1);
Pm=zeros(nt,m);
k=0;
delta=t(2)-t(1);
if round(delta)==delta
    isdaily=1;  %Interval >1day.
else
    isdaily=0;  %Invertal <24hr.
end

if isdaily
    dp=p';
    dt=t';
else
    while (i<nt)
        while ( and( round(t(i))-t(i)~=0 , i<nt))
            i=i+1;
        end
        if(i>nt)
            e=i-1;
        else
            e=i;
        end
        k=k+1;
        dt(k)=k;
        if (e-s)<=0
            disp('error');
        end
        dp(k)=sum(p(s:e))./(e-s);
        Pm(k,:)=sum(P(:,s:e),2)./(e-s);
        s=e+1;    
        i=i+1;
    end
    dp=dp(1:k);
    dt=dt(1:k);
    Pm=Pm(1:k,:);

end

end

