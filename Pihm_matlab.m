% Created by Lele Shu

%Should there be options for hourly,monthly, yearly?
%need a version control here ie v2.2 etc
%Need to add more warnings and indications of positive results

%Important user variables to change

project_name = 'vcat';
forcing_start_date = 2000;
forcing_end_date = 2009;

append_dat = true;
show_figures = 1;  %Shows figure
save_figures = 2;  %Doesnt show figure, just saves to directory
global_figures = show_figures;

%============================================
%Important USER Folders
pihm_input_dir = 'D:\\Projects\\PIHM_Matlab_Cleaned\\pihm_inputs';
pihm_output_dir = 'D:\\Projects\\PIHM_Matlab_Cleaned\\pihm_outputs';
matlab_output = 'D:\\Projects\\PIHM_Matlab_Cleaned\\matlab_output';

%============================================
disp('========================================');
disp('Only applicable to daily output');
disp('========================================');
%============================================
if (~(exist(pihm_input_dir,'dir')))
    disp('Invalid PIHM Input folder');
    return;
end
if (~(exist(pihm_output_dir,'dir')))
    disp('Invalid PIHM Output folder');
    return;
end
if (~(exist(matlab_output,'dir')))
    disp('Invalid MATLAB Output folder');
    return;
end

%============================================
%Files needed

mesh_file = strcat(pihm_input_dir, '\\',project_name,'.mesh');
att_file = strcat(pihm_input_dir, '\\',project_name,'.att');
forc_file = strcat(pihm_input_dir, '\\',project_name,'.forc');
init_file = strcat(pihm_input_dir, '\\',project_name,'.init');
soil_file = strcat(pihm_input_dir, '\\',project_name,'.soil');
geol_file = strcat(pihm_input_dir, '\\',project_name,'.geol');
riv_file = strcat(pihm_input_dir, '\\',project_name,'.riv');

if( append_dat )
    po_is_file = strcat(pihm_output_dir, '\\',project_name,'.is.dat');       %M*N; meters/day
    po_snow_file = strcat(pihm_output_dir, '\\',project_name,'.snow.dat');   %M*N; meters/day
    po_surf_file = strcat(pihm_output_dir, '\\',project_name,'.surf.dat');   %M*N; meters/day
    po_unsat_file = strcat(pihm_output_dir, '\\',project_name,'.unsat.dat'); %M*N; meters/day
    po_gw_file = strcat(pihm_output_dir, '\\',project_name,'.gw.dat');       %M*N; meters/day
    po_et0_file = strcat(pihm_output_dir, '\\',project_name,'.et0.dat');     %M*N; meters/day
    po_et1_file = strcat(pihm_output_dir, '\\',project_name,'.et1.dat');     %M*N; meters/day
    po_et2_file = strcat(pihm_output_dir, '\\',project_name,'.et2.dat');     %M*N; meters/day
    po_infil_file = strcat(pihm_output_dir, '\\',project_name,'.infil.dat'); %M*N; meters/day
    po_rech_file = strcat(pihm_output_dir, '\\',project_name,'.Rech.dat');   %M*N; meters/day
    po_qsurf_file = strcat(pihm_output_dir, '\\',project_name,'.fluxsurf.dat');
    po_qsub_file = strcat(pihm_output_dir, '\\',project_name,'.fluxsub.dat');
    po_stage_file = strcat(pihm_output_dir, '\\',project_name,'.stage.dat');
    po_rbed_file = strcat(pihm_output_dir, '\\',project_name,'.rbed.dat');
    po_rivFlx1_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx1.dat');
else
    po_is_file = strcat(pihm_output_dir, '\\',project_name,'.is');       %M*N; meters/day
    po_snow_file = strcat(pihm_output_dir, '\\',project_name,'.snow');   %M*N; meters/day
    po_surf_file = strcat(pihm_output_dir, '\\',project_name,'.surf');   %M*N; meters/day
    po_unsat_file = strcat(pihm_output_dir, '\\',project_name,'.unsat'); %M*N; meters/day
    po_gw_file = strcat(pihm_output_dir, '\\',project_name,'.gw');       %M*N; meters/day
    po_et0_file = strcat(pihm_output_dir, '\\',project_name,'.et0');     %M*N; meters/day
    po_et1_file = strcat(pihm_output_dir, '\\',project_name,'.et1');     %M*N; meters/day
    po_et2_file = strcat(pihm_output_dir, '\\',project_name,'.et2');     %M*N; meters/day
    po_infil_file = strcat(pihm_output_dir, '\\',project_name,'.infil'); %M*N; meters/day
    po_rech_file = strcat(pihm_output_dir, '\\',project_name,'.Rech');   %M*N; meters/day
    po_qsurf_file = strcat(pihm_output_dir, '\\',project_name,'.fluxsurf');
    po_qsub_file = strcat(pihm_output_dir, '\\',project_name,'.fluxsub');
    po_stage_file = strcat(pihm_output_dir, '\\',project_name,'.stage');
    po_rbed_file = strcat(pihm_output_dir, '\\',project_name,'.rbed');
    po_rivFlx1_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx1');
    
end

%============================================
%Note this is using Lele's new file formats. A modification of PIHM v2.2
names={'prcp',...   %prcp 1
    'et0', 'et1', 'et2', 'infil', 'Rech', 'fluxsurf', 'fluxsub', ... %flux 2-8
    'is', 'snow', 'surf', 'unsat', 'GW', ...    %storage 9-13
    'fluxriv1', 'fluxriv9'};    %flux of streams;   14-15

disp('Reading PIHM mesh file');

%Are the units correct?
iarea=read_area(mesh_file, global_figures, matlab_output,project_name);    %Individual area of grid.

%http://www.mathworks.com/help/matlab/ref/sum.html?searchHighlight=sum
area=sum(iarea);        %total area.

if ( area <= 0.0 )
    disp('Invalid Area generated');
    return;
end

disp('Reading PIHM Input Att file');
att = load(att_file);
[minit,rinit] = read_init(init_file,mesh_file);

disp('Reading PIHM Output files');

IS = load(po_is_file);
SNOW = load(po_snow_file);
SURF = load(po_surf_file);
UNSAT = load(po_unsat_file);
GW = load(po_gw_file);
ET0 = load(po_et0_file);
ET1 = load(po_et1_file);
ET2 = load(po_et2_file);
INFIL = load(po_infil_file);
RECH = load(po_rech_file);
QSURF = load(po_qsurf_file);
QSUB = load(po_qsub_file);
STAGE = load(po_stage_file);
RBED = load(po_rbed_file);

%============================================

disp('Creating Water Balance file');
N=length(iarea);
ts=1;
te=min([size(ET0,1),size(ET1,1),size(ET2,1),size(RECH,1),size(IS,1)]);

%Each Grid

%ET
iet0=sum(ET0(ts:te,2:N+1));
iet1=sum(ET1(ts:te,2:N+1));
iet2=sum(ET2(ts:te,2:N+1));
%Rech
ireach=sum(RECH(ts:te,2:N+1));

%Precip v2.2 forcing file .forc
% lele read precip file
% I guess this is from the new forcing file format
% How to convert from the "old" version? is the tool already here?
% Is this daily?
[p] = read_precip(forc_file);

P=p;
P=P(ts:te,:);
npsites=size(P,2);
ipm=zeros(te-ts+1,N);
patt=att(:,10)';
for i=1:npsites
    id=find(patt==i);
    ipm(:,id)=repmat(P(:,i),1,length(id));
end
ip=sum(ipm);

%Q_SURF
iqsurf=sum(QSURF(ts:te,2:N+1));
%Q_SUB
iqsub=sum(QSUB(ts:te,2:N+1));

%Load Initial Values
if ts==1
    ISinit      =	minit(:,1)';
    SNOWinit	=	minit(:,2)';
    SURFinit	=	minit(:,3)';
    GWinit      =	minit(:,5)';
    UNSATinit	=	minit(:,4)';
    STAGEinit   =   rinit(:,1)';
    RBEDinit    =   rinit(:,2)';
else
    ISinit      =	IS(ts,2:N+1);
    SNOWinit	=	SNOW(ts,2:N+1);
    SURFinit	=	SURF(ts,2:N+1);
    GWinit      =	GW(ts,2:N+1);
    UNSATinit	=	UNSAT(ts,2:N+1);
    STAGEinit   =   STAGE(ts,2:end);
    RBEDinit    =   RBED(ts,2:end);
end

%============================================
%IS
iis=IS(te,2:N+1)-ISinit;
%SNOW
isnow=SNOW(te,2:N+1)-SNOWinit;
%SURF
isurf=SURF(te,2:N+1)-SURFinit;
%UNSAT
soil=read_soil(soil_file);
isporo = soil(att(:,4),3)';
iunsat = (UNSAT(te,2:N+1)-UNSATinit).*isporo;
%GW
geol = read_geol(geol_file);
igporo=geol(att(:,4),4)';
sd=soil_depth(mesh_file)';

ids=find(GW(te,2:N+1) > sd);  % If the GW is above surface.
GW(te,ids)=sd(ids);

ids=find(GWinit>sd);  % If the GW is above surface.
GWinit(ids)=sd(ids);
igw=(GW(te,2:N+1)-GWinit).*igporo;


%============================================
if global_figures == 1
    subplot(2,1,1);
    idstg=(iis+isnow+isurf+iunsat+igw);%./iarea;
    plot(idstg);
    title('\Delta S of each grid');
    subplot(2,1,2);
    idpq=ip-(iqsurf+iqsub+iet0+iet1+iet2)./iarea;
    plot(idpq);
    
    title('P-Q-ET of each grid');
    figure(1)
    mmaps(idstg,'Delta S',ts,te, mesh_file, global_figures, matlab_output,project_name);
    figure(2)
    mmaps(idpq,'P-Q-ET',ts,te, mesh_file, global_figures, matlab_output,project_name);
    figure(3)
    mmaps(iqsurf,'Q surf accumulated',ts,te, mesh_file, global_figures, matlab_output,project_name);
    figure(4)
    mmaps(iqsub,'Q_gw',ts,te, mesh_file, global_figures, matlab_output,project_name);
    figure(5)
    mmaps(sum(RECH(ts:te,2:N+1))/area,'Recharge',ts,te, mesh_file,global_figures,  matlab_output,project_name);
elseif global_figures == 2
    subplot(2,1,1);
    idstg=(iis+isnow+isurf+iunsat+igw);%./iarea;
    plot(idstg);
    title('\Delta S of each grid');
    subplot(2,1,2);
    idpq=ip-(iqsurf+iqsub+iet0+iet1+iet2)./iarea;
    plot(idpq);
    fname = strcat(matlab_output,'\\',project_name,'_StorageByCell.png');
    print('-dpng',fname);
    
    mmaps(idstg,'Delta S',ts,te, mesh_file, global_figures, matlab_output,project_name);
    mmaps(idpq,'P-Q-ET',ts,te, mesh_file, global_figures, matlab_output,project_name);
    mmaps(iqsurf,'Q surf accumulated',ts,te, mesh_file, global_figures, matlab_output,project_name);
    mmaps(iqsub,'Q_gw',ts,te, mesh_file, global_figures, matlab_output,project_name);
    mmaps(sum(RECH(ts:te,2:N+1))/area,'Recharge',ts,te, mesh_file,global_figures,  matlab_output,project_name);
end

%============================================
fprintf('\nBalance based on each grids.\n');

IDS=(sum((iis+isnow+isurf+iunsat+igw).*iarea))/area;
TP=sum(sum(ipm).*iarea)/area;      %m
ITQ=sum((iqsurf+iqsub)./iarea);  %m3 to m;
ITE0=sum(iet0.*iarea)/area;
ITE1=sum(iet1.*iarea)/area;
ITE2=sum(iet2.*iarea)/area;

ITE=ITE0+ITE1+ITE2;

%fprintf('%s\n\n',pwd);
fprintf('DS=(sum((iis+isnow+isurf+iunsat+igw).*iarea)?/area\n');
fprintf('DS=%f\n',IDS);
fprintf('P=%f\n',TP);
fprintf('Q=%f\t \t %f\n',ITQ,ITQ/TP);
fprintf('ET=%f\t \t %f \n',ITE,ITE/TP)
fprintf('ET0=%f(%.3f)\tET1=%f(%.3f)\tET2=%f(%.3f)\n',ITE0,ITE0/TP,ITE1,ITE1/TP,ITE2,ITE2/TP);
fprintf('P-Q-ET=%f\n',TP-ITQ-ITE);
fprintf('Whole watershed with Rivers\n\n');

[Qsubwatershed,q,tq]=read_Q(po_rivFlx1_file, riv_file);

TQ=sum(sum(Qsubwatershed(ts:te,:)))/area;      %m
TP=sum(sum(ipm).*iarea)/area;      %m
TE0=sum(iet0.*iarea)/area;
TE1=sum(iet1.*iarea)/area;
TE2=sum(iet2.*iarea)/area;

TE=TE0+TE1+TE2;

%from go file start
unit = tq(2)-tq(1);
t = tq/unit/(1440/unit);
Qv2 = Qsubwatershed/area;
prcp = p(1:floor(t(end)));
%============================================
if global_figures == 1
    
    hf=figure;
    clf();
    %http://www.mathworks.com/help/matlab/ref/plotyy.html
    [ax, h1, h2] = plotyy([1:floor(t(end))],[prcp],t,Qv2,'bar','line');%,'Linewidth','10');
    
    xlabel(ax(1), 'Time (days) ');
    ylabel(ax(1), 'Precipitation (m/day) ');
    set(ax(1), 'Ydir', 'reverse');
    set(ax(1),'YAxisLocation','right')
    set(h1,'EdgeColor',[0,0,1],'FaceColor',[0,0,1]);
    set(ax(1),'ylim',[0,max([prcp;0.5])]*1.5);
    
    y=ylabel(ax(2), 'Q/Discharge (m/day)');
    set(ax(2),'YAxisLocation','left')
    set(h2, 'LineWidth', 1)
    
    ttl='Precipitation vs Discharge';
    title(ttl);
    fname = strcat(matlab_output,'\\',project_name,'_Flow.png');
    %print('-dpng',fname);
    
elseif global_figures == 2
    
    hf=figure;
    clf();
    [ax, h1, h2] = plotyy([1:floor(t(end))],[prcp],t,Qv2,'bar','line');%,'Linewidth','10');
    
    xlabel(ax(1), 'Time (days) ');
    ylabel(ax(1), 'Precipitation (m/day) ');
    set(ax(1), 'Ydir', 'reverse');
    set(ax(1),'YAxisLocation','right')
    set(h1,'EdgeColor',[0,0,1],'FaceColor',[0,0,1]);
    set(ax(1),'ylim',[0,max([prcp;0.5])]*1.5);
    
    y=ylabel(ax(2), 'Q/Discharge (m/day)');
    set(ax(2),'YAxisLocation','left')
    set(h2, 'LineWidth', 1)
    
    ttl='Precipitation vs Discharge';
    title(ttl);
    fname = strcat(matlab_output,'\\',project_name,'_Flow.png');
    print('-dpng',fname);
end

%============================================   
%from go file end

[riv,~,shp,material]=read_riv(riv_file);
shpid=riv(:,7)';
ira=zeros(1,size(riv,1));
for i=1:size(shp,1)
    tmp=shp(i,2)*shp(i,4);
    ira(shpid==i)=tmp;
end

istage=STAGE(te,2:end)-STAGEinit;
irbed=RBED(te,2:end)-RBEDinit;
istageV=sum(istage.*ira);
irbedV=sum(irbed.*ira);

%============================================
%Write file to text file
fid=fopen(strcat(matlab_output,'\\',project_name,'_iWaterbalance_',num2str(ts),'-',num2str(te),'.txt'),'w');
fprintf(fid,'Water Balance from %d, to %d\n\n',ts,te);
fprintf(fid,'\nBalance based on each grids.\n');
fprintf(fid,'\nBalance of whole watershed\n');
DS=(sum((iis+isnow+isurf+iunsat+igw).*iarea)+ istageV +irbedV)/area;

fprintf(fid,'DS=(sum((iis+isnow+isurf+iunsat+igw).*iarea)+ istageV +irbedV)/area\n');
fprintf(fid,'DS=%f\n',DS);
fprintf(fid,'P=%f\n',TP);
fprintf(fid,'Q=%f\t \tQ/P= %f\n',TQ,TQ/TP);
fprintf(fid,'ET=%f\t \tET/P= %f\n',TE,TE/TP);
fprintf(fid,'ET0=%f(%.3f)\tET1=%f(%.3f)\tET2=%f(%.3f)\n',TE0,TE0/TP,TE1,TE1/TP,TE2,TE2/TP);
fprintf(fid,'P-Q-ET=%f\n',TP-TQ-TE);

fclose(fid);
%============================================

disp('Done...');
