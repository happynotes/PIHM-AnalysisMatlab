function mmaps(data,name,ts,te,mesh_file_and_path, plot,matlab_output, project_name)

if ( nargin<1 )
    fprintf('Draw the maps of var\n');
    fprintf('Usage:\n \t mmaps(data,name,ts,te)\n\n');
    fprintf('\n\n');
    return
end
if nargin<2
    name='UNKNOWN';
end

[m,~]=size(data);

gwMax=max(max(data));
gwMin=min(min(data));
fprintf('Range of %s:\t [%f\t%f]\n',name,gwMin,gwMax);

[msh,pt]=read_mesh(mesh_file_and_path);

tri=msh(:,2:4);
x=pt(:,2);
y=pt(:,3);
z=pt(:,4);
mm=round(m/24);

if nargin==3
    if ts>=3 && ts<20
        colormap(jet(round(ts)))
    end
end
clf;
h=trisurf(tri,x,y,z,(data));

hb=colorbar;
%caxis([gwMin,(gwMax-gwMin)./2]);
%set(gca,'PlotBoxAspectRatio',[.5,.5,1])
view(2);
%set(hb, 'ylim', [0 5])

if (plot == 1)
    
    if nargin >3
        title([name,'  ',num2str(ts),'-',num2str(te)]);
        fname=strcat(matlab_output,name,num2str(ts),'-',num2str(te),'.png');
    else
        title(name);
        fname=strcat(matlab_output,name,'.png');
    end
    
    %print('-dpng',fname);
end
if (plot == 2)
    if nargin >3
        %figure('visible','off');
        title([name,'  ',num2str(ts),'-',num2str(te)]);
        fname=strcat(matlab_output,'\\',project_name,'_',name,num2str(ts),'-',num2str(te),'.png');
        %saveas(gcf, fname);
        print('-dpng',fname);
    else
        %figure('visible','off');
        title(name);
        fname=strcat(matlab_output,'\\',project_name,'_',name,'.png');
        %saveas(gcf, fname);
        print('-dpng',fname);
    end
end



hold off;
end
