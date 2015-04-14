function [ iarea ] = read_area( mesh_file_and_path, plot, save_folder,project_name )


if ( nargin ~= 4 )
    fprintf('\n[ iarea ] = read_area ( mesh_file_and_path, plot )');
    fprintf('\n\t iarea=individual area\n');
    return
end

%Lele read mesh tool
%tri is mesh returned
[tri,pt]=read_mesh(mesh_file_and_path);

x=[pt(tri(:,2),2),pt(tri(:,3),2),pt(tri(:,4),2)]';
y=[pt(tri(:,2),3),pt(tri(:,3),3),pt(tri(:,4),3)]';

%http://www.mathworks.com/help/matlab/ref/polyarea.html
iarea = polyarea(x,y);

%I thought mesh sizes were square meters
if (plot == 1)
    figure;
    hist(iarea);
    xlabel('Area km^2');
    ylabel('N');
    %print(gcf,'-dpng','histarea');
end
if (plot == 2)
    figure('visible','off');
    hist(iarea);
    xlabel('Area km^2');
    ylabel('N');
    %print(gcf,'-dpng','histarea');
    saveas(gcf, strcat(save_folder,'\\',project_name,'_areahist.png'));
end

end

