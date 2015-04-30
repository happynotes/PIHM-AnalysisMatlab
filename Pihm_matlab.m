% Created by Lele Shu

%Should there be options for hourly,monthly, yearly?
%need a version control here ie v2.2 etc
%Need to add more warnings and indications of positive results

clear;
%============================================
%Important user variables to change

project_name = 'vcat';

start_date = 2008;
end_date = 2009;

usgs_gage_filename = '01576500.txt';

append_dat = true;
append_time_stamps = true;
show_figures = 1;  %Shows figure
save_figures = 2;  %Doesnt show figure, just saves to directory
global_figures = show_figures;

plot_nashsutcliffe_e = 1; %true = 1, false = 0
%============================================
%Important USER Folders
input_output_same_location = true;
defined_input_and_output_dir='D:\\Projects\\PIHM_Matlab_Cleaned\\';
defined_input_dir='D:\\Projects\\PIHM_Matlab_Cleaned\\Input';
defined_output_dir='D:\\Projects\\PIHM_Matlab_Cleaned\\Output';

if( input_output_same_location )
    pihm_input_dir = strcat(defined_input_and_output_dir,'pihm_inputs');
    pihm_output_dir = strcat(defined_input_and_output_dir,'pihm_outputs');
    matlab_output = strcat(defined_input_and_output_dir,'matlab_output');
    usgs_input_dir = strcat(defined_input_and_output_dir,'usgs_data');
else
    pihm_input_dir = strcat(defined_input_dir,'pihm_inputs');
    pihm_output_dir = strcat(defined_output_dir,'pihm_outputs');
    matlab_output = strcat(defined_output_dir,'matlab_output');
    usgs_input_dir = strcat(defined_input_dir,'usgs_data');    
end

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
if (~(exist(usgs_input_dir,'dir')))
    disp('Invalid USGS Input folder');
    return;
end

%============================================

%============================================
%Input Files needed

mesh_file = strcat(pihm_input_dir, '\\',project_name,'.mesh');
att_file = strcat(pihm_input_dir, '\\',project_name,'.att');
forc_file = strcat(pihm_input_dir, '\\',project_name,'.forc');
init_file = strcat(pihm_input_dir, '\\',project_name,'.init');
soil_file = strcat(pihm_input_dir, '\\',project_name,'.soil');
geol_file = strcat(pihm_input_dir, '\\',project_name,'.geol');
riv_file = strcat(pihm_input_dir, '\\',project_name,'.riv');
discharge_cubic_feet_per_sec = strcat(usgs_input_dir, '\\',usgs_gage_filename);

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
    po_rivFlx0_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx0.dat');
    po_rivFlx1_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx1.dat');
    po_rivFlx2_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx2.dat');
    po_rivFlx3_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx3.dat');
    po_rivFlx4_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx4.dat');
    po_rivFlx5_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx5.dat');
    po_rivFlx6_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx6.dat');
    po_rivFlx7_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx7.dat');
    po_rivFlx8_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx8.dat');
    po_rivFlx9_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx9.dat');
   
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
    po_rivFlx0_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx0');
    po_rivFlx1_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx1');
    po_rivFlx2_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx2');
    po_rivFlx3_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx3');
    po_rivFlx4_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx4');
    po_rivFlx5_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx5');
    po_rivFlx6_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx6');
    po_rivFlx7_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx7');
    po_rivFlx8_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx8');
    po_rivFlx9_file = strcat(pihm_output_dir, '\\',project_name,'.rivFlx9');
   

end

%============================================

list_pihm_output_exts = {'is','snow','surf','unsat','gw','et0','et1','et2','infil','Rech','fluxsurf','fluxsub','stage','rbed','rivFlx0','rivFlx1'...
    'rivFlx2','rivFlx3','rivFlx4','rivFlx5','rivFlx6','rivFlx7','rivFlx8','rivFlx9'};

list_pihm_outputs = {po_is_file,po_snow_file,po_surf_file,po_unsat_file,po_gw_file,po_et0_file,po_et1_file,po_et2_file,po_infil_file,...
    po_rech_file,po_qsurf_file,po_qsub_file,po_stage_file,po_rbed_file,po_rivFlx0_file, po_rivFlx1_file,po_rivFlx2_file,po_rivFlx3_file,po_rivFlx4_file,...
    po_rivFlx5_file,po_rivFlx6_file,po_rivFlx7_file,po_rivFlx8_file,po_rivFlx9_file};
list_pihm_output_titles = {'Interception Storage [m]',...
    'Snow accumulation (water equivalent) [m]',...
    'Surface Water [m]',...
    'Soil Moisture [m]',...
    'Groundwater Head [m],'...
    'Evaporation from Canopy [m/day]',...
    'Transpiration [m/day]',...
    'Evaporation from Ground Surface [m/day]',...
    'Infiltration Rate [m/day]',...
    'Recharge  Rate [m/day]',...
    'Discharge Surface [m]',...
    'Discharge Subsurface [m]',...
    'River Stage [m],'...
    'River Bed Storage [m],'...
    'Lateral influx to the stream reach [m/day],'...
    'Lateral outflow  to the stream reach [m/day],'...
    'Surface flow to stream reach from Left terrain [m/day],'...
    'Surface flow to stream reach from right terrain [m/day],'...
    'Baseflow to stream reach from aquifer on the left [m/day],'...
    'Baseflow to stream reach from aquifer on the right [m/day],'...
    'Leakage between stream reach and river bed [m/day],'...
    'Subsurface flow to the river bed from left aquifer [m/day],'...
    'Subsurface flow to the river bed from right aquifer [m/day],'... 
    'Lateral outflux to the bed beneath river [m/day],'};
%     'Lateral influx to the bed beneath river [m/day]'};

  
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
    mmaps(idstg,'Delta S',ts,te, mesh_file, global_figures, matlab_output,project_name,append_time_stamps);
    figure(2)
    mmaps(idpq,'P-Q-ET',ts,te, mesh_file, global_figures, matlab_output,project_name,append_time_stamps);
    figure(3)
    mmaps(iqsurf,'Q surf accumulated',ts,te, mesh_file, global_figures, matlab_output,project_name,append_time_stamps);
    figure(4)
    mmaps(iqsub,'Q_gw',ts,te, mesh_file, global_figures, matlab_output,project_name,append_time_stamps);
    figure(5)
    mmaps(sum(RECH(ts:te,2:N+1))/area,'Recharge',ts,te, mesh_file,global_figures,  matlab_output,project_name,append_time_stamps);
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
    
    mmaps(idstg,'Delta S',ts,te, mesh_file, global_figures, matlab_output,project_name,append_time_stamps);
    mmaps(idpq,'P-Q-ET',ts,te, mesh_file, global_figures, matlab_output,project_name,append_time_stamps);
    mmaps(iqsurf,'Q surf accumulated',ts,te, mesh_file, global_figures, matlab_output,project_name,append_time_stamps);
    mmaps(iqsub,'Q_gw',ts,te, mesh_file, global_figures, matlab_output,project_name,append_time_stamps);
    mmaps(sum(RECH(ts:te,2:N+1))/area,'Recharge',ts,te, mesh_file,global_figures,  matlab_output,project_name,append_time_stamps);
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
    [ax, h1, h2] = plotyy((1:floor(t(end))),prcp,t,Qv2,'bar','line');%,'Linewidth','10');
    %    [ax, h1, h2] = plotyy([1:t(end)],[Qv2],t,prcp,'line','line');%,'Linewidth','10');
    xlabel(ax(1), 'Time (days) ');
    ylabel(ax(1), 'Precipitation (m/day) ');
    set(ax(1), 'Ydir', 'reverse');
    set(ax(1),'YAxisLocation','right')
    set(h1,'EdgeColor',[0,0,1],'FaceColor',[0,0,1]);
    set(ax(1),'ylim',[0,max([prcp;0.5])]*1.5);
    %    set(ax(1),'ylim',[0,max(prcp)]*0.5);
    %set(ax(1),'ylim',[0,max([prcp;0.5])]*10.5);
    %     set(ax(1),'ylim',[0,max([prcp;0.5])]);
    
    y=ylabel(ax(2), 'Q/Discharge (m/day)');
    set(ax(2),'YAxisLocation','left')
    set(h2, 'LineWidth', 1)
    set(ax(2),'ylim',[0,max(max(Qv2))*1.5]);
    
    ttl='Precipitation vs Discharge';
    if append_time_stamps == true
        formatOut = 'yyyy-mm-dd-HH-MM-SS';
        append = datestr(clock,formatOut);
        ttl = strcat( 'Precipitation vs Discharge ', append); 
    end
    
    title(ttl);
    %fname = strcat(matlab_output,'\\',project_name,'_Flow.png');
    %print('-dpng',fname);
    
elseif global_figures == 2
    
    hf=figure;
    clf();
    [ax, h1, h2] = plotyy([1:floor(t(end))],prcp,t,Qv2,'bar','line');%,'Linewidth','10');
    
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
    fname = strcat(matlab_output,'\\',project_name,'_Flow.png');
     
    if append_time_stamps == true
        formatOut = 'yyyy-mm-dd-HH-MM-SS';
        append = datestr(clock,formatOut);
        ttl = stcat( 'Precipitation vs Discharge ', append); 
        fname = strcat(matlab_output,'\\',project_name,'_Flow_',append,'.png');
    end
    
    title(ttl);
   
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
%my_wb_file_name = strcat(matlab_output,'\\',project_name,'_iWaterbalance_',num2str(ts),'-',num2str(te),'.txt');
my_wb_file_name = strcat(matlab_output,'\\',project_name,'_iWaterbalance_',append,'.txt');
if append_time_stamps == true
    formatOut = 'yyyy-mm-dd-HH-MM-SS';
    append = datestr(clock,formatOut);
    my_wb_file_name = strcat(matlab_output,'\\',project_name,'_iWaterbalance_',append,'.txt');
end

fid=fopen(my_wb_file_name,'w');
fprintf(fid,'Water Balance from %d, to %d\n\n',ts,te);
fprintf(fid,'\nBalance based on each grids.\n');
IDS=(sum((iis+isnow+isurf+iunsat+igw).*iarea))/area;
TP=sum(sum(ipm).*iarea)/area;      %m
ITQ=sum((iqsurf+iqsub)./iarea);  %m3 to m;
ITE0=sum(iet0.*iarea)/area;
ITE1=sum(iet1.*iarea)/area;
ITE2=sum(iet2.*iarea)/area;
ITE=ITE0+ITE1+ITE2;

fprintf(fid,'DS=(sum((iis+isnow+isurf+iunsat+igw).*iarea)?/area\n');
fprintf(fid,'DS=%f\n',IDS);
fprintf(fid,'P=%f\n',TP);
fprintf(fid,'Q=%f\t \t %f\n',ITQ,ITQ/TP);
fprintf(fid,'ET=%f\t \t %f \n',ITE,ITE/TP);
fprintf(fid,'ET0=%f(%.3f)\tET1=%f(%.3f)\tET2=%f(%.3f)\n',ITE0,ITE0/TP,ITE1,ITE1/TP,ITE2,ITE2/TP);
fprintf(fid,'P-Q-ET=%f\n',TP-ITQ-ITE);
fprintf(fid,'Whole watershed with Rivers\n\n');

fprintf(fid,'***************************************************************************************.\n');
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
A = size(list_pihm_output_titles);
len = A(2);
for n = 1:len
    clear title;
    clear variable_input_filename_path;
    clear precipitation_input_filename_path;
    clear ext;
    clear output_filename_path;

    
    title = cell2mat(list_pihm_output_titles(n));
    variable_input_filename_path = cell2mat(list_pihm_outputs(n));
    precipitation_input_filename_path = forc_file;
    ext = cell2mat(list_pihm_output_exts(n));
    output_filename_path = strcat(matlab_output, '\\',project_name,'_model_',ext,'.png');
         
    if append_time_stamps == true
        formatOut = 'yyyy-mm-dd-HH-MM-SS';
        append = datestr(clock,formatOut);
        title = strcat(title,' ', append);
        output_filename_path = strcat(matlab_output, '\\',project_name,'_model_',append,'_', ext,'.png');
    end
    
    plot_pihm_variables(variable_input_filename_path, precipitation_input_filename_path,output_filename_path,title);
end


%============================================
if exist(discharge_cubic_feet_per_sec, 'file')
    %Load USGS observed data, NOTE DATA MUST BE DAILY
    fileID = fopen(discharge_cubic_feet_per_sec,'r');
    formatSpec = '%*s %*s %*s %d %*s';
    allflow = fscanf(fileID,formatSpec);
    fclose(fileID);
    
    Q=Qsubwatershed/area;
    
    ufactor=24*3600*.348^3;
    obs = allflow(t,1)*ufactor/area;   %convert to m/day from cms;
    
    %http://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient
    %GOAL IS 1
    [NSE,id,MatchedData]=nashsutcliffe([t,obs],[t,Q]);
    
    %============================================
    hf=figure;
    clf();
    % vl=yline(yline<t(end));
    % if ~isempty(vl)
    %     for i=1:length(vl)
    %         line([vl(i),vl(i)],[0,1],'LineWidth',4,'Color',[.8 .8 .8]);hold on;
    %     end
    % end
    
    if ( plot_nashsutcliffe_e)
        subplot(4,1,[1,3]);
        [ax, h1, h2] = plotyy([1:t(end)],[Q,obs],t,prcp,'line','bar');%,'Linewidth','10');
        xlabel(ax(2), 'Time (days) ');
        ylabel(ax(2), 'Precipitation (m/day) ');
        set(ax(2), 'Ydir', 'reverse');
        set(ax(2),'YAxisLocation','right')
        set(h2,'FaceColor',[0.5,0.5,0.5],'EdgeColor',[0,0,1]);
        ymax=ceil(max(prcp)*1.5*10)/10;
        set(ax(2),'ylim',[0,ymax]);
        set(ax(2),'ytick',0:ymax/5:ymax);
        
        ymax=ceil(max(max(Q),max(obs))*1.1*100)/100;
        set(ax(1),'YAxisLocation','left')
        set(ax(1),'ylim',[0,ymax]);
        set(ax(1),'ytick',0:ymax/5:ymax);
        ylabel(ax(1), 'Q/Discharge (m/day)');
        set(h1(1),'Color', 'r');    %Simulation
        set(h1(2),'Color', 'b');    %Observation
        annotation('textbox', [0.1, 0.40, 0, 0], 'string', start_date);
        annotation('textbox', [0.9, 0.40, 0, 0], 'string', end_date);
        nse_title = strcat('NSE =  ',num2str(NSE));
        annotation('textbox', [0.75, 0.40, 0, 0], 'string', nse_title);
        clear title local_ttl
        local_ttl='Precipitation vs Discharge';
        
        %Output Files
        usgs_vs_model_output = strcat(matlab_output, '\\',project_name,'_model_vs_usgs_discharge.png');
        
        if append_time_stamps == true
            formatOut = 'yyyy-mm-dd-HH-MM-SS';
            append = datestr(clock,formatOut);
            local_ttl = strcat(local_ttl,' ', append);
            usgs_vs_model_output = strcat(matlab_output, '\\',project_name,'_model_vs_usgs_discharge',append,'.png');
        end
        
        title(local_ttl);
        
        leg=legend([h2,h1(1),h1(2)],'Precipitation','Simulation','Observation','location','southoutside');
        set(leg,'FontSize',10);
        
        subplot(4,1,4);
        plot(t,MatchedData(:,2) - MatchedData(:,3),'color','r');
        local_ttl='Nashsutcliffe-E';
        annotation('textbox', [0.45, 0.05, 0.0, 0.0], 'string', local_ttl);
        
        % hold on;
        % plot(t,MatchedData(:,3),'color','b');
        if( global_figures == 2 )
            print('-dpng',usgs_vs_model_output);
        end
        
        hold off;
    else
        [ax, h1, h2] = plotyy([1:t(end)],[Q,obs],t,prcp,'line','bar');%,'Linewidth','10');
        xlabel(ax(2), 'Time (days) ');
        ylabel(ax(2), 'Precipitation (m/day) ');
        set(ax(2), 'Ydir', 'reverse');
        set(ax(2),'YAxisLocation','right')
        set(h2,'FaceColor',[0.5,0.5,0.5],'EdgeColor',[0,0,1]);
        ymax=ceil(max(prcp)*1.5*10)/10;
        set(ax(2),'ylim',[0,ymax]);
        set(ax(2),'ytick',0:ymax/5:ymax);
        
        ymax=ceil(max(max(Q),max(obs))*1.1*100)/100;
        set(ax(1),'YAxisLocation','left')
        set(ax(1),'ylim',[0,ymax]);
        set(ax(1),'ytick',0:ymax/5:ymax);
        ylabel(ax(1), 'Q/Discharge (m/day)');
        set(h1(1),'Color', 'r');    %Simulation
        set(h1(2),'Color', 'b');    %Observation
        annotation('textbox', [0.1, 0.20, 0, 0], 'string', start_date);
        annotation('textbox', [0.9, 0.20, 0, 0], 'string', end_date);
        nse_title = strcat('NSE =  ',num2str(NSE));
        annotation('textbox', [0.75, 0.20, 0, 0], 'string', nse_title);
        clear title local_ttl
        local_ttl='Precipitation vs Discharge';
        
        usgs_vs_model_output = strcat(matlab_output, '\\',project_name,'_model_vs_usgs_discharge.png');
        
        if append_time_stamps == true
            formatOut = 'yyyy-mm-dd-HH-MM-SS';
            append = datestr(clock,formatOut);
            local_ttl = strcat(local_ttl,' ', append);
            usgs_vs_model_output = strcat(matlab_output, '\\',project_name,'_model_vs_usgs_discharge',append,'.png');
        end
        
        title(local_ttl);
        leg=legend([h2,h1(1),h1(2)],'Precipitation','Simulation','Observation','location','southoutside');
        set(leg,'FontSize',10);
        
        if( global_figures == 2 )
            print('-dpng',usgs_vs_model_output);
        end
        
        hold off;
    end
    
    
    fprintf('NSE=\t%g\n\n\n',NSE);
    
end
%============================================


disp('Done...');
