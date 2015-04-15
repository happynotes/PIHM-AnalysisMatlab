function [] = plot_pihm_variables( variable_input_filename_path, precipitation_input_filename_path, output_filename_path, my_title )

if ~exist(variable_input_filename_path, 'file')
    display( strcat('Missing file: ', variable_input_filename_path ));
    return
end

[mm] = load(variable_input_filename_path);
[m,n]=size(mm);
ts=1;
te=m;

figure;

subplot(2,1,1);
title(my_title);
line(ts:te,mm(ts:te,2:end));
set(gca,'XTick',[]);


subplot(2,1,2);
tmean=mean(mm(ts:te,2:end),2);
h=plot(ts:te,tmean);
set(h,'color','r','linewidth',1);
ylabel('Average')
allmean=mean(tmean);
hold on;
hl=line([ts,te],[allmean,allmean]);
set(hl,'color',[0,0,1]);

hold on
t=mm(:,1)/1440;

[ p,dt,Pm ] = read_precip(precipitation_input_filename_path);

pend=t(end);
tp=ts:te;
prcp=p(ts:te);
if (max(tmean)>min(tmean))
    set(gca,'ylim',[min(tmean),max(tmean)*1.5]);
end
ax1 = gca;
hold on

set(gca,'yDir','reverse');
legend('Mean');

x = ts:te;
plotyy(x,allmean,x,prcp)

%hold on
%ylabel(ax(1), 'Precipitation (m/day) ');

%ax2 = axes('Position',get(ax1,'Position'),'YAxisLocation','right','Color','none','XColor','k','YColor','k');
%
% linkaxes([ax1 ax2],'x');
% %%plotprcp(prcp,tp);
% 
% T=1:length(prcp);
% 
% if sum(size(prcp)-size(T))~=0
%     fprintf('ERROR:\n\tsize of Prcp: [%d,%d]\n',size(prcp));
%     fprintf('\tsize of Time: [%d,%d]\n',size(T));
%     fprintf('\tStop ploting\n\n');
%     return;
% end
% 
% 
% hold on;
% h1 = gca;
% h2 = axes('Position',get(h1,'Position'));
% 
% bar(T,prcp,'r','EdgeColor','r')
% line(T,prcp,'Color','b','LineWidth',4);
% set(h2,'YAxisLocation','right','yDir','reverse','Color','none','XTickLabel',[])
% set(h2,'XLim',get(h1,'XLim'),'Layer','bottom')
% 
% ts=T(1);
% te=T(end);
% hl=line([ts,te],[mean(prcp),mean(prcp)]);
% set(hl,'color',[0,0,1]);
% set(gca,'ylim',[0,max(prcp)*1.5]);
% 
% set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
% ylabel('Precipitation (m/day)');


print('-dpng',output_filename_path);
hold off;

end

