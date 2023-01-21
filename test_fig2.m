% Script to plot Figure 2
% BY Yangkang Chen
% 
% Initialized: Jan, 2022
% Revised:     Jan, 2023
% This script takes about 1-2 minutes
% 
%% Please first download the MATseisdl package
% svn co https://github.com/chenyk1990/MATseisdl/trunk MATseisdl

clc;clear;close all;

addpath(genpath('./MATseisdl'))
addpath(genpath('./subroutines'))

%
%% point 
x0=32*0.01;
y0=24*0.01;
z0=36*0.01;

%% line
y1=linspace(0.3,0.7,40);
x1=0.3*ones(size(y1));
z1=0.15*ones(size(x1));

[x2,y2]=meshgrid([0:0.1:1],[0:0.1:1]);

z2=0*ones(size(x2));

z2=[0:11]*0.03;
x2=0*ones(size(z2));
y2=0*ones(size(z2));

z22=[0:30]*0.03;
x22=0*ones(size(z22))+0.05;
y22=0*ones(size(z22));

figure('units','normalized','Position',[0.2 0.4 0.4, 0.5],'color','w');
plot3(y0,x0,z0,'rp','linewidth',1,'MarkerSize',16,'MarkerFaceColor','r');hold on;
xlabel('Y (km)','Fontsize',16,'Fontweight','normal');
ylabel('X (km)','Fontsize',16,'Fontweight','normal');
zlabel('Depth (km)','Fontsize',16,'Fontweight','normal');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','normal');
view(gca,[-59.7 29.3991240875912]);
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
grid on

% plot3(y1,x1,z1,'r-','linewidth',3);
% plot3(y1,x1,z1,'ro','linewidth',3,'MarkerSize',12);
xlim([0,1]);ylim([0,1]);zlim([0,1]);

plot3(y2,x2,z2,'bv','linewidth',3,'MarkerSize',12);
plot3(y22,x22,z22,'gv','linewidth',3,'MarkerSize',12);
% xlim([0,1]);ylim([0,1]);zlim([0,0.2]);

set(get(gca,'xlabel'),'rotation',35,'VerticalAlignment','middle')
set(get(gca,'ylabel'),'rotation',-10,'VerticalAlignment','top')

%% add text
text(0,0.6,{'Vp=3.5 km/s','Vs=2.0 km/s','\rho=1000 g/m^3','fm=30 Hz','Source (sx,sy,sz)=(0.32,0.24,0.36) km','Grids (nx,ny,nz)=(100,100,100)','Sampling (dx,dy,dz)=(0.01,0.01,0.01) km'},'color','k','Fontsize',16,'fontweight','normal');


a1=axes('Parent',gcf,'Position',[0.1,0.8,0.2,0.2]);
mm=[
 0,1,0,
 1,0,0,
 0,0,0];
M=[0,0,0,1,0,0];
Mw = 5.0;
focalmech(M,0,0,Mw,'k','xyz');
axis equal
axis off
% scatter(x,y);box on;set(gca,'xticklabel',[]);ylabel('North (m)');set(gca,'Linewidth',1.5,'Fontsize',12);

% Create arrow
annotation(gcf,'arrow',[0.347 0.25],[0.565 0.811],'linewidth',2);


print(gcf,'-depsc','-r300','fig2.eps');
print(gcf,'-dpng','-r300','fig2.png');

