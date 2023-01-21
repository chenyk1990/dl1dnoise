% Script to plot Figure 6
% BY Yangkang Chen
% 
% Initialized: Jan, 2022
% Revised:     Jan, 2023
% This script takes about 1-2 minutes
% 
%% Please first download the MATseisdl package
% svn co https://github.com/chenyk1990/MATseisdl/trunk MATseisdl

clc;clear;close all;

addpath(genpath('./MATseisdl'));
addpath(genpath('./subroutines'));


load('data/syn1.mat');
% original script and datapath
% rvz=zeros(351,12);
% rvx=zeros(351,12);
% rvy=zeros(351,12);
% rsf_read(rvz,'/Users/chenyk/chenyk/micro_dl2/syn3d/rvz2.rsf');
% rsf_read(rvx,'/Users/chenyk/chenyk/micro_dl2/syn3d/rvx2.rsf');
% rsf_read(rvy,'/Users/chenyk/chenyk/micro_dl2/syn3d/rvy2.rsf');

rvz=rvz*1000000000;
rvx=rvx*1000000000;
rvy=rvy*1000000000;


t=[0:350]*0.001+0.15;
x=[1:12];
clip=0.2;
dc=[rvz,rvx,rvy];
randn('state',20202122);
dn=dc+0.02*randn(size(dc));
%%

%% patch size l1*l2
l1=32;l2=1;

c1=l1;c2=64;%size of the 1D cosine dictionary (if c2>c1, overcomplete)
%% DCT dictionary (dctmtx will generates orthogonal transform)
dct=zeros(c1,c2);
for k=0:1:c2-1
    V=cos([0:1:c1-1]'*k*pi/c2);
    if k>0
        V=V-mean(V);
    end
    dct(:,k+1)=V/norm(V);
end
D=dct;
% DCT=kron(dct,dct);%2D DCT dictionary (64,256)

%% decompose the image into patches:
X=dl_patch(dn,1,l1,1,l1/2,1);


%% OMP using DCT
nd=size(X,2);
K=3;
ph=1;
tic
for i2=1:nd
    G(:,i2)=dl_omp0(D,X(:,i2),K);
end
toc

%further constrain it to be sparser
G=dl_pthresh(G,'ph',ph);
X2=D*G;

[n1,n2]=size(dn);
d2=dl_patch_inv(X2,1,n1,n2,l1,1,l1/2,1);

%% KSVD
param.T=K;      %sparsity level
param.D=D;    %initial D
param.niter=30; %number of K-SVD iterations to perform; default: 10
param.mode=1;   %1: sparsity; 0: error
param.K=c2;     %number of atoms, dictionary size
tic
[Dksvd,Gksvd]=dl_ksvd(X,param);
toc
Gksvd0=Gksvd;
Gksvd=dl_pthresh(Gksvd0,'ph',ph);
X1=Dksvd*Gksvd;
d1=dl_patch_inv(X1,1,n1,n2,l1,1,l1/2,1);

%% SGK
param.T=K;      %sparsity level
param.D=D;    %initial D
param.niter=30; %number of K-SVD iterations to perform; default: 10
param.mode=1;   %1: sparsity; 0: error
param.K=c2;     %number of atoms, dictionary size
tic
[Dsgk,Gsgk]=dl_sgk(X,param);
toc
Gsgk0=Gsgk;
Gsgk=dl_pthresh(Gsgk0,'ph',ph);
X11=Dsgk*Gsgk;
d11=dl_patch_inv(X11,1,n1,n2,l1,1,l1/2,1);

% figure('units','normalized','Position',[0.2 0.4 1 0.8]);
% for ia=1:64
%     subplot(8,8,ia);plot([1:l1],D(:,ia),'k','linewidth',2);
% 
%     if ia==1
%         ylim([-0.5,0.5]);
%     end
%     set(gca,'Linewidth',2.0,'Fontsize',10);
%     ytickformat('%.1f');
%     if ismember(ia,1:8:64)
%         ylabel('Amplitude','Fontsize',10);
%     else
%         set(gca,'yticklabel',[]);
% 
%     end
% 
%     if ismember(ia,57:64)
%         xlabel('Sample NO','Fontsize',10);
%     else
%         set(gca,'xticklabel',[]);
%     end
%     axis off;
%     xlim([1,l1]);
% end
% annotation(gcf,'rectangle',[0.126 0.101 0.782 0.834],'linewidth',2);
% print(gcf,'-depsc','-r300','syn1_atom0_new.eps');

% figure('units','normalized','Position',[0.2 0.4 1 0.8]);
% for ia=1:64
%     subplot(8,8,ia);plot([1:l1],Dksvd(:,ia),'k','linewidth',2);
% 
%     if ia==1
%         ylim([-0.5,0.5]);
%     end
%     if ismember(ia,[1,2,3,4,5,6,7,8,9,10,11,12,33])
%         subplot(8,8,ia);plot([1:l1],Dksvd(:,ia),'r','linewidth',2);
%     end
%     set(gca,'Linewidth',2.0,'Fontsize',10);
%     ytickformat('%.1f');
%     if ismember(ia,1:8:64)
%         ylabel('Amplitude','Fontsize',10);
%     else
%         set(gca,'yticklabel',[]);
% 
%     end
% 
%     if ismember(ia,57:64)
%         xlabel('Sample NO','Fontsize',10);
%     else
%         set(gca,'xticklabel',[]);
%     end
%     axis off;
%     xlim([1,l1]);
% end
% annotation(gcf,'rectangle',[0.126 0.101 0.782 0.834],'linewidth',2);
% print(gcf,'-depsc','-r300','syn1_atom1_new.eps');

% figure('units','normalized','Position',[0.2 0.4 1 0.8]);
% for ia=1:64
%     subplot(8,8,ia);plot([1:l1],Dsgk(:,ia),'k','linewidth',2);
% 
%     if ia==1
%         ylim([-0.5,0.5]);
%     end
%     if ismember(ia,[4,5,6,7,10,11,30,44,59])
%         subplot(8,8,ia);plot([1:l1],Dsgk(:,ia),'r','linewidth',2);
%     end
%     set(gca,'Linewidth',2.0,'Fontsize',10);
%     ytickformat('%.1f');
%     if ismember(ia,1:8:64)
%         ylabel('Amplitude','Fontsize',10);
%     else
%         set(gca,'yticklabel',[]);
% 
%     end
% 
%     if ismember(ia,57:64)
%         xlabel('Sample NO','Fontsize',10);
%     else
%         set(gca,'xticklabel',[]);
%     end
%     axis off
%     xlim([1,l1]);
% end
% annotation(gcf,'rectangle',[0.126 0.101 0.782 0.834],'linewidth',2);
% print(gcf,'-depsc','-r300','syn1_atom2_new.eps');

D1=D;
Dksvd1=Dksvd;
Dsgk1=Dsgk;


%%% Second case
load('data/syn2.mat');
% original script and datapath
% rvz=zeros(351,12);
% rvx=zeros(351,12);
% rvy=zeros(351,12);
% rsf_read(rvz,'/Users/chenyk/chenyk/micro_dl2/syn3d/rvz2.rsf');
% rsf_read(rvx,'/Users/chenyk/chenyk/micro_dl2/syn3d/rvx2.rsf');
% rsf_read(rvy,'/Users/chenyk/chenyk/micro_dl2/syn3d/rvy2.rsf');

rvz=rvz*1000000000;
rvx=rvx*1000000000;
rvy=rvy*1000000000;

%
% rvz=zeros(501,100);
% rvx=zeros(501,100);
% rvy=zeros(501,100);
%
% rsf_read(rvz,'../syn3d/rvz.rsf');rvz=rvz*1000000000;
% rsf_read(rvx,'../syn3d/rvx.rsf');rvx=rvx*1000000000;
% rsf_read(rvy,'../syn3d/rvy.rsf');rvy=rvy*1000000000;

[n1,n2]=size(rvz);
t=[0:n1-1]*0.001;
x=[1:n2];
clip=0.2;
dc=[rvz,rvx,rvy];
randn('state',20202122);
dn=dc+0.02*randn(size(dc));
%%

%% patch size l1*l2
l1=32;l2=1;
l1=64;
c1=l1;c2=64;
%% DCT dictionary (dctmtx will generates orthogonal transform)
dct=zeros(c1,c2);
for k=0:1:c2-1
    V=cos([0:1:c1-1]'*k*pi/c2);
    if k>0
        V=V-mean(V);
    end
    dct(:,k+1)=V/norm(V);
end
D=dct;
% DCT=kron(dct,dct);%2D DCT dictionary (64,256)

%% plot the first 64 atoms
% figure;
% for ia=1:16
%     subplot(4,4,ia);plot(dct(:,ia));
% end

%% decompose the image into patches:
X=dl_patch(dn,1,l1,1,l1/2,1);


%% OMP using DCT
nd=size(X,2);
K=3;%sparsity level
ph=1;
tic
for i2=1:nd
    G(:,i2)=dl_omp0(D,X(:,i2),K);
end
toc

%further constrain it to be sparser
G=dl_pthresh(G,'ph',ph);
X2=D*G;

[n1,n2]=size(dn);
d2=dl_patch_inv(X2,1,n1,n2,l1,1,l1/2,1);

%% KSVD
param.T=K;      %sparsity level
param.D=D;    %initial D
param.niter=30; %number of K-SVD iterations to perform; default: 10
param.mode=1;   %1: sparsity; 0: error
param.K=c2;     %number of atoms, dictionary size
tic
[Dksvd,Gksvd]=dl_ksvd(X,param);
toc
Gksvd0=Gksvd;
Gksvd=dl_pthresh(Gksvd0,'ph',ph);
X1=Dksvd*Gksvd;
d1=dl_patch_inv(X1,1,n1,n2,l1,1,l1/2,1);

%% SGK
param.T=K;      %sparsity level
param.D=D;    %initial D
param.niter=30; %number of K-SVD iterations to perform; default: 10
param.mode=1;   %1: sparsity; 0: error
param.K=c2;     %number of atoms, dictionary size
tic
[Dsgk,Gsgk]=dl_sgk(X,param);
toc
Gsgk0=Gsgk;
Gsgk=dl_pthresh(Gsgk0,'ph',ph);
X11=Dsgk*Gsgk;
d11=dl_patch_inv(X11,1,n1,n2,l1,1,l1/2,1);

% figure('units','normalized','Position',[0.2 0.4 1 0.8],'color','w');
% for ia=1:64
%     subplot(8,8,ia);plot([1:l1],D(:,ia),'k','linewidth',2);
% 
%     if ia==1
%         ylim([-0.5,0.5]);
%     end
%     set(gca,'Linewidth',2.0,'Fontsize',10);
%     ytickformat('%.1f');
%     if ismember(ia,1:8:64)
%         ylabel('Amplitude','Fontsize',10);
%     else
%         set(gca,'yticklabel',[]);
% 
%     end
% 
%     if ismember(ia,57:64)
%         xlabel('Sample NO','Fontsize',10);
%     else
%         set(gca,'xticklabel',[]);
%     end
%     axis off;
%     xlim([1,l1]);
% end
% annotation(gcf,'rectangle',[0.126 0.101 0.782 0.834],'linewidth',2);
% print(gcf,'-depsc','-r300','syn2_atom0_new.eps');

% figure('units','normalized','Position',[0.2 0.4 1 0.8],'color','w');
% for ia=1:64
%     subplot(8,8,ia);plot([1:l1],Dksvd(:,ia),'k','linewidth',2);
% 
%     if ia==1
%         ylim([-0.5,0.5]);
%     end
%     if ismember(ia,[1,2,3,4,5,6,7,8,9,10,12,13,15,16,22,25,38,44,47,51])
%         subplot(8,8,ia);plot([1:l1],Dksvd(:,ia),'r','linewidth',2);
%     end
%     set(gca,'Linewidth',2.0,'Fontsize',10);
%     ytickformat('%.1f');
%     if ismember(ia,1:8:64)
%         ylabel('Amplitude','Fontsize',10);
%     else
%         set(gca,'yticklabel',[]);
% 
%     end
% 
%     if ismember(ia,57:64)
%         xlabel('Sample NO','Fontsize',10);
%     else
%         set(gca,'xticklabel',[]);
%     end
%     axis off;
%     xlim([1,l1]);
% end
% annotation(gcf,'rectangle',[0.126 0.101 0.782 0.834],'linewidth',2);
% print(gcf,'-depsc','-r300','syn2_atom1_new.eps');

% figure('units','normalized','Position',[0.2 0.4 1 0.8],'color','w');
% for ia=1:64
%     subplot(8,8,ia);plot([1:l1],Dsgk(:,ia),'k','linewidth',2);
% 
%     if ia==1
%         ylim([-0.5,0.5]);
%     end
%     if ismember(ia,[2,3,4,5,6,7,8,9,10,46,58])
%         subplot(8,8,ia);plot([1:l1],Dsgk(:,ia),'r','linewidth',2);
%     end
%     set(gca,'Linewidth',2.0,'Fontsize',10);
%     ytickformat('%.1f');
%     if ismember(ia,1:8:64)
%         ylabel('Amplitude','Fontsize',10);
%     else
%         set(gca,'yticklabel',[]);
% 
%     end
% 
%     if ismember(ia,57:64)
%         xlabel('Sample NO','Fontsize',10);
%     else
%         set(gca,'xticklabel',[]);
%     end
%     axis off;
%     xlim([1,l1]);
% end
% annotation(gcf,'rectangle',[0.126 0.101 0.782 0.834],'linewidth',2);
% print(gcf,'-depsc','-r300','syn2_atom2_new.eps');

%% combined

% save fig6.mat
% load fig6.mat

figure('units','normalized','Position',[0.0 0.0 0.45, 1],'color','w');

x0=0.05;y0=0.05;dy=0;dx=0;
%length: 0.5x0.5, 0.5x0.25, 0.25x0.5
nr=24;%number of stations in the first column
nc=16;
xm=0.05; %%middle blank
%% axis XY
dh=(1-y0*2)/nr;%dw=0.37;
dw=(1-x0*2-xm)/nc;
dh1=dh;%axis height
t=[0:1999]*0.0005;

%% P(1,1)
l1=32;
for ir=1:8
    for ic=1:8
        a1=axes('Parent',gcf,'Position',[x0+(ic-1)*dw,y0+(ir-1)*dh,dw,dh]);
        ia=(8-ir)*8+ic;
        if ismember(ia,[4,5,6,7,10,11,30,44,59])
            plot([1:l1],Dsgk1(:,ia),'r','linewidth',2);
        else
            plot([1:l1],Dsgk1(:,ia),'k','linewidth',2);
        end
        xlim([1-3,l1+3]);ylim([-0.5-0.1,0.5+0.1]);
        axis off;
    end
end

%% P(1,2)
for ir=1:8
    for ic=1:8
        a1=axes('Parent',gcf,'Position',[x0+(ic-1)*dw,y0+(ir+8-1)*dh,dw,dh]);
        ia=(8-ir)*8+ic;
        if ismember(ia,[1,2,3,4,5,6,7,8,9,10,11,12,33])
            plot([1:l1],Dksvd1(:,ia),'r','linewidth',2);
        else
            plot([1:l1],Dksvd1(:,ia),'k','linewidth',2);
        end

        xlim([1-3,l1+3]);ylim([-0.5-0.1,0.5+0.1]);
        axis off;
    end
end

%% P(1,3)
for ir=1:8
    for ic=1:8
        a1=axes('Parent',gcf,'Position',[x0+(ic-1)*dw,y0+(ir+8*2-1)*dh,dw,dh]);
        ia=(8-ir)*8+ic;
        plot([1:l1],D1(:,ia),'k','linewidth',2);
        xlim([1-3,l1+3]);ylim([-0.5-0.1,0.5+0.1]);
        axis off;
    end
end


l1=64;
%% P(2,1)
for ir=1:8
    for ic=1:8
        a1=axes('Parent',gcf,'Position',[x0+(ic+8-1)*dw+xm,y0+(ir-1)*dh,dw,dh]);
        ia=(8-ir)*8+ic;
        if ismember(ia,[2,3,4,5,6,7,8,9,10,46,58])
            plot([1:l1],Dsgk(:,ia),'r','linewidth',2);
        else
            plot([1:l1],Dsgk(:,ia),'k','linewidth',2);
        end
        xlim([1-3,l1+3]);ylim([-0.5-0.1,0.5+0.1]);
        axis off;
    end
end

%% P(2,2)
for ir=1:8
    for ic=1:8
        a1=axes('Parent',gcf,'Position',[x0+(ic+8-1)*dw+xm,y0+(ir+8-1)*dh,dw,dh]);
        ia=(8-ir)*8+ic;
        if ismember(ia,[1,2,3,4,5,6,7,8,9,10,12,13,15,16,22,25,38,44,47,51])
            plot([1:l1],Dksvd(:,ia),'r','linewidth',2);
        else
            plot([1:l1],Dksvd(:,ia),'k','linewidth',2);
        end

        xlim([1-3,l1+3]);ylim([-0.5-0.1,0.5+0.1]);
        axis off;
    end
end

%% P(2,3)
for ir=1:8
    for ic=1:8
        a1=axes('Parent',gcf,'Position',[x0+(ic+8-1)*dw+xm,y0+(ir+8*2-1)*dh,dw,dh]);
        ia=(8-ir)*8+ic;
        plot([1:l1],D(:,ia),'k','linewidth',2);
        xlim([1-3,l1+3]);ylim([-0.5-0.1,0.5+0.1]);
        axis off;
    end
end

annotation(gcf,'rectangle',[x0 y0 8*dw 24*dh],'linewidth',2,'color','k');
annotation(gcf,'rectangle',[x0+8*dw+xm y0 8*dw 24*dh],'linewidth',2,'color','k');


annotation(gcf,'rectangle',[x0 y0 8*dw 16*dh],'linewidth',2,'color','k');
annotation(gcf,'rectangle',[x0 y0 8*dw 8*dh],'linewidth',2,'color','k');
annotation(gcf,'rectangle',[x0+8*dw+xm y0 8*dw 16*dh],'linewidth',2,'color','k');
annotation(gcf,'rectangle',[x0+8*dw+xm y0 8*dw 8*dh],'linewidth',2,'color','k');
% annotation(gcf,'rectangle',[x0 y0 1-x0*2 1-y0*2],'linewidth',2);

annotation(gcf,'textbox',[0 y0/2+24*dh 0.05 0.05],'String','(a)','Fontsize',20,'fontweight','bold','FitBoxToText','off','EdgeColor','none');
annotation(gcf,'textbox',[0 y0/2+16*dh 0.05 0.05],'String','(c)','Fontsize',20,'fontweight','bold','FitBoxToText','off','EdgeColor','none');
annotation(gcf,'textbox',[0 y0/2+8*dh 0.05 0.05],'String','(e)','Fontsize',20,'fontweight','bold','FitBoxToText','off','EdgeColor','none');

annotation(gcf,'textbox',[x0+8*dw y0/2+24*dh 0.05 0.05],'String','(b)','Fontsize',20,'fontweight','bold','FitBoxToText','off','EdgeColor','none');
annotation(gcf,'textbox',[x0+8*dw y0/2+16*dh 0.05 0.05],'String','(d)','Fontsize',20,'fontweight','bold','FitBoxToText','off','EdgeColor','none');
annotation(gcf,'textbox',[x0+8*dw y0/2+8*dh 0.05 0.05],'String','(f)','Fontsize',20,'fontweight','bold','FitBoxToText','off','EdgeColor','none');


print(gcf,'-depsc','-r300','fig6.eps');
print(gcf,'-dpng','-r300','fig6.png');







