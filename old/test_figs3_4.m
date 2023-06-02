% Script to plot Figures 3 and 4
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
% figure('units','normalized','Position',[0.2 0.4 0.5, 0.4],'color','w');
% subplot(1,6,1:2);dl_wigb(rvz,1,x,t,clip);xlabel('Channel');title('Vz');ylabel('Time (s)');
% subplot(1,6,3:4);dl_wigb(rvx,1,x,t,clip);xlabel('Channel');title('Vx');
% subplot(1,6,5:6);dl_wigb(rvy,1,x,t,clip);xlabel('Channel');title('Vy');

dc=[rvz,rvx,rvy];
randn('state',20202122);
dn=dc+0.02*randn(size(dc));
% figure;dl_imagesc([dc,dn]);
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

%% plot the first 64 atoms
% figure;
% for ia=1:16
%     subplot(4,4,ia);plot(dct(:,ia));
% end

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
% figure;dl_imagesc([dn,d2,dn-d2]);

% figure('units','normalized');
% imagesc(G);colormap(jet);colorbar;caxis([-0.5,0.5]);%colorbar;
% ylabel('Atom NO','Fontsize',16);
% xlabel('Patch NO','Fontsize',16);
% title('Coefficients Matrix','Fontsize',16);
% set(gca,'Linewidth',1.5,'Fontsize',16);
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
% figure;dl_imagesc([dn,d1,dn-d1]);

% figure('units','normalized','Position',[0.2 0.4 0.8, 0.8],'color','w');
% for ia=1:64
%     subplot(8,8,ia);plot(Dksvd(:,ia));
% end
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
% figure;dl_imagesc([dn,d11,dn-d11]);
dl_snr(dc,d11)

% figure('units','normalized','Position',[0.2 0.4 1 0.8]);
% for ia=1:64
%     subplot(8,8,ia);plot([1:l1],D(:,ia),'b','linewidth',2);
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
%     
%     xlim([1,l1]);
% end
% % print(gcf,'-depsc','-r300','syn1_atom0.eps');
% 
% figure('units','normalized','Position',[0.2 0.4 1 0.8]);
% for ia=1:64
%     subplot(8,8,ia);plot([1:l1],Dksvd(:,ia),'b','linewidth',2);
%     
%     if ia==1
%         ylim([-0.5,0.5]);
%     end
%         if ismember(ia,[1,2,3,4,5,6,7,8,9,10,11,12,33])
%            subplot(8,8,ia);plot([1:l1],Dksvd(:,ia),'r','linewidth',2);
%         end
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
%     
%     xlim([1,l1]);
% end
% print(gcf,'-depsc','-r300','syn1_atom1.eps');

% figure('units','normalized','Position',[0.2 0.4 1 0.8]);
% for ia=1:64
%     subplot(8,8,ia);plot([1:l1],Dsgk(:,ia),'b','linewidth',2);
%     
%     if ia==1
%         ylim([-0.5,0.5]);
%     end
%         if ismember(ia,[4,5,6,7,10,11,30,44,59])
%            subplot(8,8,ia);plot([1:l1],Dsgk(:,ia),'r','linewidth',2);
%         end
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
%     
%     xlim([1,l1]);
% end
% print(gcf,'-depsc','-r300','syn1_atom2.eps');



% figure('units','normalized','Position',[0.2 0.4 0.5, 0.8],'color','w');
% subplot(2,6,1:2);dl_wigb(rvz,1,x,t,clip);title('Vz');ylabel('Time (s)');%xlabel('Channel');
% subplot(2,6,3:4);dl_wigb(rvx,1,x,t,clip);title('Vx');set(gca,'ytick',[]);%xlabel('Channel');
% subplot(2,6,5:6);dl_wigb(rvy,1,x,t,clip);title('Vy');set(gca,'ytick',[]);%xlabel('Channel');
% 
% subplot(2,6,7:8);dl_wigb(dn(:,1:12),1,x,t,clip);xlabel('Channel');title('Vz');ylabel('Time (s)');
% subplot(2,6,9:10);dl_wigb(dn(:,13:24),1,x,t,clip);xlabel('Channel');title('Vx');set(gca,'ytick',[])
% subplot(2,6,11:12);dl_wigb(dn(:,25:36),1,x,t,clip);xlabel('Channel');title('Vy');set(gca,'ytick',[])
% print(gcf,'-depsc','-r300','fig2.eps');

tic
d3=dl_bp(dn,0.004,0,0.1,40,45,0);
toc

inds=50:151;%because the SNR is not always the best metric
dl_snr(dc(inds,:),dn(inds,:))
dl_snr(dc(inds,:),d1(inds,:))
dl_snr(dc(inds,:),d11(inds,:))
dl_snr(dc(inds,:),d2(inds,:))
dl_snr(dc(inds,:),d3(inds,:))


%% for zooming
x1=8.5;x2=11.5;y1=0.2;y2=0.3;

figure('units','normalized','Position',[0.2 0.4 0.5, 1],'color','w');
subplot(3,4,1:2);dl_wigb([rvz,rvx,rvy],1,1:36,t(:),clip);%title('Vz');ylabel('Time (s)');%xlabel('Channel');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
ylabel('Time (s)','Fontsize',12,'fontweight','normal');
title('Clean data','Fontsize',12,'fontweight','normal');
text(-5,0.1,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(6,0.16,'Vz','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(18,0.16,'Vx','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(30,0.16,'Vy','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
fp1=dl_ap2fp([10,0.2]);fp2=dl_ap2fp([11,0.24]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([14,0.33]);fp2=dl_ap2fp([13,0.3]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
dl_framebox(x1,x2,y1,y2,'g',2);

subplot(3,4,3:4);dl_wigb(dn,1,1:36,t(:),clip);%title('Vz');ylabel('Time (s)');%xlabel('Channel');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
title('Noisy data','Fontsize',12,'fontweight','normal');
text(-5,0.1,'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(6,0.16,'Vz','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(18,0.16,'Vx','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(30,0.16,'Vy','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
fp1=dl_ap2fp([10,0.2]);fp2=dl_ap2fp([11,0.24]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([14,0.33]);fp2=dl_ap2fp([13,0.3]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
dl_framebox(x1,x2,y1,y2,'g',2);

subplot(3,4,5:6);dl_wigb(d3,1,1:36,t(:),clip);%title('Vz');ylabel('Time (s)');%xlabel('Channel');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
ylabel('Time (s)','Fontsize',12,'fontweight','normal');
title('Denoised (BP)','Fontsize',12,'fontweight','normal');
text(-5,0.1,'(c)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(6,0.16,'Vz','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(18,0.16,'Vx','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(30,0.16,'Vy','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
fp1=dl_ap2fp([10,0.2]);fp2=dl_ap2fp([11,0.24]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([14,0.33]);fp2=dl_ap2fp([13,0.3]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
dl_framebox(x1,x2,y1,y2,'g',2);

subplot(3,4,7:8);dl_wigb(d2,1,1:36,t(:),clip);%title('Vz');ylabel('Time (s)');%xlabel('Channel');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
title('Denoised (DCT)','Fontsize',12,'fontweight','normal');
text(-5,0.1,'(d)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(6,0.16,'Vz','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(18,0.16,'Vx','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(30,0.16,'Vy','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
fp1=dl_ap2fp([10,0.2]);fp2=dl_ap2fp([11,0.24]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([14,0.33]);fp2=dl_ap2fp([13,0.3]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
dl_framebox(x1,x2,y1,y2,'g',2);

subplot(3,4,9:10);dl_wigb(d1,1,1:36,t(:),clip);%title('Vz');ylabel('Time (s)');%xlabel('Channel');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
ylabel('Time (s)','Fontsize',12,'fontweight','normal');
xlabel('Channel','Fontsize',12,'fontweight','normal');
title('Denoised (KSVD)','Fontsize',12,'fontweight','normal');
text(-5,0.1,'(e)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(6,0.16,'Vz','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(18,0.16,'Vx','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(30,0.16,'Vy','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
fp1=dl_ap2fp([10,0.2]);fp2=dl_ap2fp([11,0.24]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([14,0.33]);fp2=dl_ap2fp([13,0.3]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
dl_framebox(x1,x2,y1,y2,'g',2);

subplot(3,4,11:12);dl_wigb(d11,1,1:36,t(:),clip);%title('Vz');ylabel('Time (s)');%xlabel('Channel');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
xlabel('Channel','Fontsize',12,'fontweight','normal');
title('Denoised (SGK)','Fontsize',12,'fontweight','normal');
text(-5,0.1,'(f)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(6,0.16,'Vz','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(18,0.16,'Vx','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(30,0.16,'Vy','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
fp1=dl_ap2fp([10,0.2]);fp2=dl_ap2fp([11,0.24]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([14,0.33]);fp2=dl_ap2fp([13,0.3]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
dl_framebox(x1,x2,y1,y2,'g',2);
print(gcf,'-depsc','-r300','fig3.eps');
print(gcf,'-dpng','-r300','fig3.png');


figure('units','normalized','Position',[0.2 0.4 0.5, 1],'color','w');
subplot(3,4,1:2);dl_wigb_two([rvz,rvx,rvy],1,1:36,t(:),clip,[rvz,rvx,rvy]);%title('Vz');ylabel('Time (s)');%xlabel('Channel');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
ylabel('Time (s)','Fontsize',12,'fontweight','normal');
title('Clean data','Fontsize',12,'fontweight','normal');
text(x1+(x2-x1)/38*(-5),y1+(y2-y1)/0.25*(0.1-0.15),'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlim([x1,x2]);ylim([y1,y2]);

subplot(3,4,3:4);dl_wigb_two(dn,1,1:36,t(:),clip,[rvz,rvx,rvy]);%title('Vz');ylabel('Time (s)');%xlabel('Channel');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
title('Noisy data','Fontsize',12,'fontweight','normal');
text(x1+(x2-x1)/38*(-5),y1+(y2-y1)/0.25*(0.1-0.15),'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlim([x1,x2]);ylim([y1,y2]);

subplot(3,4,5:6);dl_wigb_two(d3,1,1:36,t(:),clip,[rvz,rvx,rvy]);%title('Vz');ylabel('Time (s)');%xlabel('Channel');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
ylabel('Time (s)','Fontsize',12,'fontweight','normal');
title('Denoised (BP)','Fontsize',12,'fontweight','normal');
text(x1+(x2-x1)/38*(-5),y1+(y2-y1)/0.25*(0.1-0.15),'(c)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlim([x1,x2]);ylim([y1,y2]);

subplot(3,4,7:8);dl_wigb_two(d2,1,1:36,t(:),clip,[rvz,rvx,rvy]);%title('Vz');ylabel('Time (s)');%xlabel('Channel');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
title('Denoised (DCT)','Fontsize',12,'fontweight','normal');
text(x1+(x2-x1)/38*(-5),y1+(y2-y1)/0.25*(0.1-0.15),'(d)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlim([x1,x2]);ylim([y1,y2]);

subplot(3,4,9:10);dl_wigb_two(d1,1,1:36,t(:),clip,[rvz,rvx,rvy]);%title('Vz');ylabel('Time (s)');%xlabel('Channel');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
ylabel('Time (s)','Fontsize',12,'fontweight','normal');
xlabel('Channel','Fontsize',12,'fontweight','normal');
title('Denoised (KSVD)','Fontsize',12,'fontweight','normal');
text(x1+(x2-x1)/38*(-5),y1+(y2-y1)/0.25*(0.1-0.15),'(e)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlim([x1,x2]);ylim([y1,y2]);

subplot(3,4,11:12);dl_wigb_two(d11,1,1:36,t(:),clip,[rvz,rvx,rvy]);%title('Vz');ylabel('Time (s)');%xlabel('Channel');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
xlabel('Channel','Fontsize',12,'fontweight','normal');
title('Denoised (SGK)','Fontsize',12,'fontweight','normal');
text(x1+(x2-x1)/38*(-5),y1+(y2-y1)/0.25*(0.1-0.15),'(f)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
xlim([x1,x2]);ylim([y1,y2]);
print(gcf,'-depsc','-r300','fig4.eps');
print(gcf,'-dpng','-r300','fig4.png');


