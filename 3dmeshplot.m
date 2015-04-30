%This is the code for 3d view of the PIHM mesh 
filename=importdata('projectName.txt');
readname=strcat(filename,'.mesh');
fid = importdata(readname{1});

Ele=fid(1,1);
Node=fid(1,2);
t=zeros(Ele,3);
p=zeros(Node,3);
nbrs=zeros(Ele,3);

cdata=zeros(Ele,1);

for i=1:Ele
    for j=1:3
        t(i,j)=fid(i+1,j+1);
        nbrs(i,j)=fid(i+1,j+4);
    end
end
for i=1:Node
    for j=1:2
        p(i,j)=fid(1+Ele+i,j+1);
    end
    p(i,3)=fid(1+Ele+i,5);
end


tr = TriRep(t,p);
nbrs = neighbors(tr);

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);

axis equal;
trisurf(tr,'EdgeColor',[0.2 0.2 0.2],'facecolor','interp');
hold on;

zlimits = [min(p(:,3)) max(p(:,3))];
demcmap(zlimits);

pp=p;
for i=1:Node
    pp(i,3)=fid(1+Ele+i,5)-500;
end
trbed=TriRep(t,pp);
trisurf(trbed,'EdgeColor','none','LineStyle','none','FaceColor',[0.5 0.5 0.6]);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca,'ZTick',[]);
vector2=zeros(2,1);
for i=1:Ele
    for k=1:3
        if(isnan(nbrs(i,k)))
            vector2(1)=mod(k,3)+1;
            vector2(2)=mod(k+1,3)+1;
        
p1 = [p(t(i,vector2(1)),1) p(t(i,vector2(1)),2) p(t(i,vector2(1)),3)];
p2 = [pp(t(i,vector2(1)),1) pp(t(i,vector2(1)),2) (pp(t(i,vector2(1)),3)+p(t(i,vector2(1)),3))/2];
p3 = [pp(t(i,vector2(2)),1) pp(t(i,vector2(2)),2) (pp(t(i,vector2(2)),3)+p(t(i,vector2(2)),3))/2];
p4 = [p(t(i,vector2(2)),1) p(t(i,vector2(2)),2) p(t(i,vector2(2)),3)]; 
p5 = [pp(t(i,vector2(1)),1) pp(t(i,vector2(1)),2) pp(t(i,vector2(1)),3)];
p6 = [pp(t(i,vector2(2)),1) pp(t(i,vector2(2)),2) pp(t(i,vector2(2)),3)];
x = [p1(1) p2(1) p3(1) p4(1)];
y = [p1(2) p2(2) p3(2) p4(2)];
z = [p1(3) p2(3) p3(3) p4(3)];

u = [p2(1) p5(1) p6(1) p3(1)];
v = [p2(2) p5(2) p6(2) p3(2)];
w = [p2(3) p5(3) p6(3) p3(3)];
fill3(x, y, z, [161 137 125]/255,'EdgeColor','none');
fill3(u,v,w,[137 153 189]/255,'EdgeColor','none');
        end
    end
end
axis off;

az=-57.5;
el=30;
view(az,el);
