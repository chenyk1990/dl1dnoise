% Script to plot Figure 5
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

load('data/syn2.mat');

% original script and datapath
% rvz=zeros(501,100);
% rvx=zeros(501,100);
% rvy=zeros(501,100);
% rsf_read(rvz,'/Users/chenyk/chenyk/micro_dl2/syn3d/rvz.rsf');
% rsf_read(rvx,'/Users/chenyk/chenyk/micro_dl2/syn3d/rvx.rsf');
% rsf_read(rvy,'/Users/chenyk/chenyk/micro_dl2/syn3d/rvy.rsf');

rvz=rvz*1000000000;
rvx=rvx*1000000000;
rvy=rvy*1000000000;

[n1,n2]=size(rvz);
t=[0:n1-1]*0.001;
x=[1:n2];
clip=0.2;
figure('units','normalized','Position',[0.2 0.4 0.5, 0.4],'color','w');
subplot(1,6,1:2);dl_wigb(rvz,1,x,t,clip);xlabel('Channel');title('Vz');ylabel('Time (s)');
subplot(1,6,3:4);dl_wigb(rvx,1,x,t,clip);xlabel('Channel');title('Vx');
subplot(1,6,5:6);dl_wigb(rvy,1,x,t,clip);xlabel('Channel');title('Vy');

dc=[rvz,rvx,rvy];
randn('state',20202122);
dn=dc+0.02*randn(size(dc));
figure;dl_imagesc([dc,dn]);
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
figure;dl_imagesc([dn,d2,dn-d2]);dl_snr(dc,d2)

figure('units','normalized');
imagesc(G);colormap(jet);colorbar;caxis([-0.5,0.5]);%colorbar;
ylabel('Atom NO','Fontsize',16);
xlabel('Patch NO','Fontsize',16);
title('Coefficients Matrix','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);
% print(gcf,'-depsc','-r300','real_G.eps');


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
figure;dl_imagesc([dn,d1,dn-d1]);dl_snr(dc,d1)

figure('units','normalized','Position',[0.2 0.4 0.8, 0.8],'color','w');
for ia=1:64
    subplot(8,8,ia);plot(Dksvd(:,ia));
end
% print(gcf,'-depsc','-r300','syn_atom_ksvd.eps');

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
figure;dl_imagesc([dn,d11,dn-d11]);dl_snr(dc,d11)

figure('units','normalized','Position',[0.2 0.4 0.8, 0.8],'color','w');
for ia=1:64
    subplot(8,8,ia);plot(Dsgk(:,ia));
end
% print(gcf,'-depsc','-r300','syn_atom_sgk.eps');

figure('units','normalized','Position',[0.2 0.4 0.8, 0.8],'color','w');
for ia=1:64
    subplot(8,8,ia);plot(D(:,ia));
end
% print(gcf,'-depsc','-r300','syn_atom_dct.eps');

dl_snr(dc,dn)
dl_snr(dc,d1)
dl_snr(dc,d11)
dl_snr(dc,d2)


figure('units','normalized','Position',[0.2 0.4 1 0.8]);
for ia=1:64
    subplot(8,8,ia);plot([1:l1],D(:,ia),'b','linewidth',2);
    
    if ia==1
        ylim([-0.5,0.5]);
    end
    set(gca,'Linewidth',2.0,'Fontsize',10);
    ytickformat('%.1f');
    if ismember(ia,1:8:64)
        ylabel('Amplitude','Fontsize',10);
    else
        set(gca,'yticklabel',[]);
        
    end
    
    if ismember(ia,57:64)
        xlabel('Sample NO','Fontsize',10);
    else
        set(gca,'xticklabel',[]);
    end
    
    xlim([1,l1]);
end
% print(gcf,'-depsc','-r300','syn2_atom0.eps');

figure('units','normalized','Position',[0.2 0.4 1 0.8]);
for ia=1:64
    subplot(8,8,ia);plot([1:l1],Dksvd(:,ia),'b','linewidth',2);
    
    if ia==1
        ylim([-0.5,0.5]);
    end
        if ismember(ia,[1,2,3,4,5,6,7,8,9,10,12,13,15,16,22,25,38,44,47,51])
           subplot(8,8,ia);plot([1:l1],Dksvd(:,ia),'r','linewidth',2);
        end
    set(gca,'Linewidth',2.0,'Fontsize',10);
    ytickformat('%.1f');
    if ismember(ia,1:8:64)
        ylabel('Amplitude','Fontsize',10);
    else
        set(gca,'yticklabel',[]);
        
    end
    
    if ismember(ia,57:64)
        xlabel('Sample NO','Fontsize',10);
    else
        set(gca,'xticklabel',[]);
    end
    
    xlim([1,l1]);
end
% print(gcf,'-depsc','-r300','syn2_atom1.eps');

figure('units','normalized','Position',[0.2 0.4 1 0.8]);
for ia=1:64
    subplot(8,8,ia);plot([1:l1],Dsgk(:,ia),'b','linewidth',2);
    
    if ia==1
        ylim([-0.5,0.5]);
    end
        if ismember(ia,[2,3,4,5,6,7,8,9,10,46,58])
           subplot(8,8,ia);plot([1:l1],Dsgk(:,ia),'r','linewidth',2);
        end
    set(gca,'Linewidth',2.0,'Fontsize',10);
    ytickformat('%.1f');
    if ismember(ia,1:8:64)
        ylabel('Amplitude','Fontsize',10);
    else
        set(gca,'yticklabel',[]);
        
    end
    
    if ismember(ia,57:64)
        xlabel('Sample NO','Fontsize',10);
    else
        set(gca,'xticklabel',[]);
    end
    
    xlim([1,l1]);
end
% print(gcf,'-depsc','-r300','syn2_atom2.eps');




figure('units','normalized','Position',[0.2 0.4 0.55, 0.8],'color','w');
subplot(3,1,1);dl_imagesc([dn,d2,dn-d2]);title('DCT, SNR=7.65 dB (0.81 s), l1=64, K=3, ph=1','Fontsize',12,'fontweight','normal');
subplot(3,1,2);dl_imagesc([dn,d1,dn-d1]);title('KSVD, SNR=12.97 dB (26.81 s), l1=64, K=3, ph=1','Fontsize',12,'fontweight','normal');
subplot(3,1,3);dl_imagesc([dn,d11,dn-d11]);title('SGK, SNR=10.22 dB (2.01 s), l1=64, K=3, ph=1','Fontsize',12,'fontweight','normal');
% print(gcf,'-depsc','-r300','fig3old.eps');


%% better figure
x=1:300;
figure('units','normalized','Position',[0.2 0.4 0.55, 1],'color','w');
subplot(3,3,1);dl_imagesc([dc],0.2,2,x,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
text(-60,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Noise-free data','Fontsize',12,'fontweight','normal');
text(50,0.03,'Vz','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(150,0.03,'Vx','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(250,0.03,'Vy','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
% text(10,0.425,'SNR=3.56 dB','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','left');

subplot(3,3,2);dl_imagesc([dn],0.2,2,x,t);set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');%ylabel('Time (s)','Fontsize',12,'fontweight','normal');
text(-60,-0.05,'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Noisy data','Fontsize',12,'fontweight','normal');
text(50,0.03,'Vz','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(150,0.03,'Vx','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(250,0.03,'Vy','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(10,0.425,'SNR=3.56 dB','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','left');

subplot(3,3,4);dl_imagesc(d2,0.2,2,x,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
text(-60,-0.05,'(c)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Denoised (DCT)','Fontsize',12,'fontweight','normal');
text(50,0.03,'Vz','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(150,0.03,'Vx','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(250,0.03,'Vy','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(10,0.45,{'SNR=7.65 dB','Cost=0.81 s'},'color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','left');
fp1=dl_ap2fp([235-200,0.06]);fp2=dl_ap2fp([245-200,0.11]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([235-100,0.06]);fp2=dl_ap2fp([245-100,0.11]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([235,0.06]);fp2=dl_ap2fp([245,0.11]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');

subplot(3,3,5);dl_imagesc(d1,0.2,2,x,t);set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
text(-60,-0.05,'(d)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Denoised (KSVD)','Fontsize',12,'fontweight','normal');
text(50,0.03,'Vz','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(150,0.03,'Vx','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(250,0.03,'Vy','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(10,0.45,{'SNR=12.97 dB','Cost=26.81 s'},'color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','left');
fp1=dl_ap2fp([235-200,0.06]);fp2=dl_ap2fp([245-200,0.11]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([235-100,0.06]);fp2=dl_ap2fp([245-100,0.11]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([235,0.06]);fp2=dl_ap2fp([245,0.11]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');

subplot(3,3,6);dl_imagesc(d11,0.2,2,x,t);set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
text(-60,-0.05,'(e)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Denoised (SGK)','Fontsize',12,'fontweight','normal');
text(50,0.03,'Vz','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(150,0.03,'Vx','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(250,0.03,'Vy','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(10,0.45,{'SNR=10.22 dB','Cost=2.01 s'},'color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','left');
fp1=dl_ap2fp([235-200,0.06]);fp2=dl_ap2fp([245-200,0.11]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([235-100,0.06]);fp2=dl_ap2fp([245-100,0.11]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([235,0.06]);fp2=dl_ap2fp([245,0.11]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');

subplot(3,3,7);dl_imagesc(dn-d2,0.2,2,x,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');xlabel('Channel','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
text(-60,-0.05,'(f)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Noise (DCT)','Fontsize',12,'fontweight','normal');
text(50,0.03,'Vz','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(150,0.03,'Vx','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(250,0.03,'Vy','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
fp1=dl_ap2fp([235-200,0.1]);fp2=dl_ap2fp([245-200,0.15]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([235-100,0.1]);fp2=dl_ap2fp([245-100,0.15]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([235,0.1]);fp2=dl_ap2fp([245,0.15]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');


subplot(3,3,8);dl_imagesc(dn-d1,0.2,2,x,t);xlabel('Channel','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
text(-60,-0.05,'(g)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Noise (KSVD)','Fontsize',12,'fontweight','normal');
text(50,0.03,'Vz','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(150,0.03,'Vx','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(250,0.03,'Vy','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
fp1=dl_ap2fp([235-200,0.1]);fp2=dl_ap2fp([245-200,0.15]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([235-100,0.1]);fp2=dl_ap2fp([245-100,0.15]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([235,0.1]);fp2=dl_ap2fp([245,0.15]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');

subplot(3,3,9);dl_imagesc(dn-d11,0.2,2,x,t);xlabel('Channel','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
text(-60,-0.05,'(h)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Noise (SGK)','Fontsize',12,'fontweight','normal');
text(50,0.03,'Vz','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(150,0.03,'Vx','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
text(250,0.03,'Vy','color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','center');
fp1=dl_ap2fp([235-200,0.1]);fp2=dl_ap2fp([245-200,0.15]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([235-100,0.1]);fp2=dl_ap2fp([245-100,0.15]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([235,0.1]);fp2=dl_ap2fp([245,0.15]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
print(gcf,'-depsc','-r300','fig5.eps');
print(gcf,'-dpng','-r300','fig5.png');


% print(gcf,'-dpng','-r300','result1.png');

% 
% figure('units','normalized','Position',[0.2 0.4 0.5, 0.8],'color','w');
% subplot(2,6,1:2);dl_wigb(rvz,1,x,t,clip);title('Vz');ylabel('Time (s)');%xlabel('Channel');
% subplot(2,6,3:4);dl_wigb(rvx,1,x,t,clip);title('Vx');set(gca,'ytick',[]);%xlabel('Channel');
% subplot(2,6,5:6);dl_wigb(rvy,1,x,t,clip);title('Vy');set(gca,'ytick',[]);%xlabel('Channel');
% 
% subplot(2,6,7:8);dl_wigb(dn(:,1:12),1,x,t,clip);xlabel('Channel');title('Vz');ylabel('Time (s)');
% subplot(2,6,9:10);dl_wigb(dn(:,13:24),1,x,t,clip);xlabel('Channel');title('Vx');set(gca,'ytick',[])
% subplot(2,6,11:12);dl_wigb(dn(:,25:36),1,x,t,clip);xlabel('Channel');title('Vy');set(gca,'ytick',[])
% print(gcf,'-depsc','-r300','fig3.eps');
% 
% figure('units','normalized','Position',[0.2 0.4 0.5, 0.8],'color','w');
% subplot(2,6,1:2);dl_wigb(d2(:,1:12),1,x,t,clip);title('Vz');ylabel('Time (s)');%xlabel('Channel');
% subplot(2,6,3:4);dl_wigb(d2(:,13:24),1,x,t,clip);title('Vx');set(gca,'ytick',[]);%xlabel('Channel');
% subplot(2,6,5:6);dl_wigb(d2(:,25:36),1,x,t,clip);title('Vy');set(gca,'ytick',[]);%xlabel('Channel');
% 
% subplot(2,6,7:8);dl_wigb(d1(:,1:12),1,x,t,clip);xlabel('Channel');title('Vz');ylabel('Time (s)');
% subplot(2,6,9:10);dl_wigb(d1(:,13:24),1,x,t,clip);xlabel('Channel');title('Vx');set(gca,'ytick',[])
% subplot(2,6,11:12);dl_wigb(d1(:,25:36),1,x,t,clip);xlabel('Channel');title('Vy');set(gca,'ytick',[])
% print(gcf,'-depsc','-r300','fig33.eps');
% 
% figure('units','normalized','Position',[0.2 0.4 0.5, 0.8],'color','w');
% subplot(2,6,1:2);dl_wigb(dn(:,1:12)-d2(:,1:12),1,x,t,clip);title('Vz');ylabel('Time (s)');%xlabel('Channel');
% subplot(2,6,3:4);dl_wigb(dn(:,13:24)-d2(:,13:24),1,x,t,clip);title('Vx');set(gca,'ytick',[]);%xlabel('Channel');
% subplot(2,6,5:6);dl_wigb(dn(:,25:36)-d2(:,25:36),1,x,t,clip);title('Vy');set(gca,'ytick',[]);%xlabel('Channel');
% 
% subplot(2,6,7:8);dl_wigb(dn(:,1:12)-d1(:,1:12),1,x,t,clip);xlabel('Channel');title('Vz');ylabel('Time (s)');
% subplot(2,6,9:10);dl_wigb(dn(:,13:24)-d1(:,13:24),1,x,t,clip);xlabel('Channel');title('Vx');set(gca,'ytick',[])
% subplot(2,6,11:12);dl_wigb(dn(:,25:36)-d1(:,25:36),1,x,t,clip);xlabel('Channel');title('Vy');set(gca,'ytick',[])
% print(gcf,'-depsc','-r300','fig333.eps');
% 
