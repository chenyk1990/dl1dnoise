% Script to plot Figures 8,9,10
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

load data/real1.mat

%%indices for a different threshold (applied on the coefficients)
inds=[134:219,257:330,456:462,661:694];


%% patch size l1*l2
l1=64;l2=1;

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
X=dl_patch(d,1,l1,1,l1/2,1);


%% OMP using DCT
nd=size(X,2);
K=3;
ph=0.5;
tic
for i2=1:nd
    G(:,i2)=dl_omp0(D,X(:,i2),K);
end
toc

%further constrain it to be sparser
G=dl_pthresh(G,'ph',ph);
X2=D*G;

[n1,n2]=size(d);
d2=dl_patch_inv(X2,1,n1,n2,l1,1,l1/2,1);
% figure;dl_imagesc([d,d2,d-d2]);

% figure('units','normalized');
% imagesc(G);colormap(jet);colorbar;caxis([-0.5,0.5]);%colorbar;
% ylabel('Atom NO','Fontsize',16);
% xlabel('Patch NO','Fontsize',16);
% title('Coefficients Matrix','Fontsize',16);
% set(gca,'Linewidth',1.5,'Fontsize',16);


%% K-SVD
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
[n1,n2]=size(d);
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
% print(gcf,'-depsc','-r300','real_atom0.eps');

% figure('units','normalized','Position',[0.2 0.4 1 0.8]);
% for ia=1:64
%     subplot(8,8,ia);plot([1:l1],Dksvd(:,ia),'b','linewidth',2);
%     
%     if ia==1
%         ylim([-0.5,0.5]);
%     end
%         if ismember(ia,[12,13,14,18,53,54,57,58])
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
% print(gcf,'-depsc','-r300','real_atom1.eps');

% figure('units','normalized','Position',[0.2 0.4 1 0.8]);
% for ia=1:64
%     subplot(8,8,ia);plot([1:l1],Dsgk(:,ia),'b','linewidth',2);
%     
%     if ia==1
%         ylim([-0.5,0.5]);
%     end
%         if ismember(ia,[6,7,8,9,10,11,12,13,14,32,47,52,54,55])
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
% print(gcf,'-depsc','-r300','real_atom2.eps');


%% a secondary QC of the results
%DCT
XX=dl_patch(d(:,inds),1,l1,1,l1/2,1);%for the special group of traces
ph2=1.0;
ndd=size(XX,2);
for i2=1:ndd
    GG(:,i2)=dl_omp0(D,XX(:,i2),K);
end
GG=dl_pthresh(GG,'ph',ph2);
XX2=D*GG;
[nn1,nn2]=size(d(:,inds));
dd2=dl_patch_inv(XX2,1,nn1,nn2,l1,1,l1/2,1);
d2(:,inds)=dd2;
%KSVD
for i2=1:ndd
    GGksvd(:,i2)=dl_omp0(Dksvd,XX(:,i2),K);
end
GGksvd=dl_pthresh(GGksvd,'ph',ph2);
XX1=Dksvd*GGksvd;
dd1=dl_patch_inv(XX1,1,nn1,nn2,l1,1,l1/2,1);
d1(:,inds)=dd1;
%SGK
for i2=1:ndd
    GGsgk(:,i2)=dl_omp0(Dsgk,XX(:,i2),K);
end
GGsgk=dl_pthresh(GGsgk,'ph',ph2);
XX11=Dsgk*GGsgk;
dd11=dl_patch_inv(XX11,1,nn1,nn2,l1,1,l1/2,1);
d11(:,inds)=dd11;
%% a secondary QC of the results


% figure;dl_imagesc([d,d1,d-d1]);
% figure;dl_imagesc([d,d11,d-d11]);
% figure;dl_imagesc([d,d2,d-d2]);

dn=d;
% figure('units','normalized','Position',[0.2 0.4 0.55, 1],'color','w');
% subplot(3,1,1);dl_imagesc([dn,d2,dn-d2]);title('DCT,  (4.48 s), l1=64, K=3, ph=0.5');
% subplot(3,1,2);dl_imagesc([dn,d1,dn-d1]);title('KSVD,  (169.71 s), l1=64, K=3, ph=0.5');
% subplot(3,1,3);dl_imagesc([dn,d11,dn-d11]);title('SGK, (14.00 s), l1=64, K=3, ph=0.5');
% print(gcf,'-depsc','-r300','fig8.eps');

% save test_fig8.mat d d1 d11 d2

% figure;
% subplot(2,1,1);plot(d(:,1));
% subplot(2,1,2);plot(d1(:,1));
% 
% ix=800;
% figure;
% subplot(2,1,1);plot(d(:,ix));
% subplot(2,1,2);plot(d1(:,ix));
% 
% ix=400;
% figure;
% subplot(3,1,1);plot(d(:,ix));
% subplot(3,1,2);plot(d1(:,ix));
% subplot(3,1,3);plot(d2(:,ix));
% 
% ix=900;
% figure;
% subplot(3,1,1);plot(d(:,ix));
% subplot(3,1,2);plot(d1(:,ix));
% subplot(3,1,3);plot(d2(:,ix));
% 
% ix=910;
% figure;
% subplot(4,1,1);plot(d(:,ix));ylim([-999,1200]);title('Noisy');
% subplot(4,1,2);plot(d1(:,ix));ylim([-999,1200]);title('KSVD');
% subplot(4,1,3);plot(d11(:,ix));ylim([-999,1200]);title('SGK');
% subplot(4,1,4);plot(d2(:,ix));ylim([-999,1200]);title('DCT');
% print(gcf,'-depsc','-r300','fig88.eps');


t=[0:1000]*0.002;
x=1:918;
figure('units','normalized','Position',[0.2 0.4 0.8, 1],'color','w');
subplot(3,3,2);dl_imagesc([dn],2000,2,x,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
text(-120,-0.15,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Raw data','Fontsize',12,'fontweight','normal');
fp1=dl_ap2fp([580,0.3]);fp2=dl_ap2fp([625,0.5]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([690,0.3]);fp2=dl_ap2fp([735,0.5]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([785,0.25]);fp2=dl_ap2fp([835,0.45]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');

subplot(3,3,4);dl_imagesc(d2,2000,2,x,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
text(-120,-0.15,'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Denoised (DCT)','Fontsize',12,'fontweight','normal');
text(50,1.75,{'Cost=4.48 s'},'color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','left');
fp1=dl_ap2fp([580,0.3]);fp2=dl_ap2fp([625,0.5]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([690,0.3]);fp2=dl_ap2fp([735,0.5]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([785,0.25]);fp2=dl_ap2fp([835,0.45]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');

subplot(3,3,5);dl_imagesc(d1,2000,2,x,t);set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
text(-120,-0.15,'(c)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Denoised (KSVD)','Fontsize',12,'fontweight','normal');
text(50,1.75,{'Cost=169.71 s'},'color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','left');
fp1=dl_ap2fp([580,0.3]);fp2=dl_ap2fp([625,0.5]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([690,0.3]);fp2=dl_ap2fp([735,0.5]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([785,0.25]);fp2=dl_ap2fp([835,0.45]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');

subplot(3,3,6);dl_imagesc(d11,2000,2,x,t);set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
text(-120,-0.15,'(d)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Denoised (SGK)','Fontsize',12,'fontweight','normal');
text(50,1.75,{'Cost=14.00 s'},'color','k','Fontsize',12,'fontweight','normal','HorizontalAlignment','left');
fp1=dl_ap2fp([580,0.3]);fp2=dl_ap2fp([625,0.5]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([690,0.3]);fp2=dl_ap2fp([735,0.5]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([785,0.25]);fp2=dl_ap2fp([835,0.45]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');

subplot(3,3,7);dl_imagesc(dn-d2,2000,2,x,t);ylabel('Time (s)','Fontsize',12,'fontweight','normal');xlabel('Channel','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
text(-120,-0.15,'(e)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Noise (DCT)','Fontsize',12,'fontweight','normal');
fp1=dl_ap2fp([40,0.3]);fp2=dl_ap2fp([52,0.57]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([560,0.40]);fp2=dl_ap2fp([580,0.65]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([785,0.25]);fp2=dl_ap2fp([850,0.55]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');

subplot(3,3,8);dl_imagesc(dn-d1,2000,2,x,t);xlabel('Channel','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
text(-120,-0.15,'(f)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Noise (KSVD)','Fontsize',12,'fontweight','normal');
fp1=dl_ap2fp([40,0.3]);fp2=dl_ap2fp([52,0.57]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([560,0.40]);fp2=dl_ap2fp([580,0.65]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([785,0.25]);fp2=dl_ap2fp([850,0.55]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');


subplot(3,3,9);dl_imagesc(dn-d11,2000,2,x,t);xlabel('Channel','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','normal');
text(-120,-0.15,'(g)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
title('Noise (SGK)','Fontsize',12,'fontweight','normal');
fp1=dl_ap2fp([40,0.3]);fp2=dl_ap2fp([52,0.57]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([560,0.40]);fp2=dl_ap2fp([580,0.65]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');
fp1=dl_ap2fp([785,0.25]);fp2=dl_ap2fp([850,0.55]);annotation('arrow',[fp1(1),fp2(1)],[fp1(2),fp2(2)],'linewidth',2,'color','r');

print(gcf,'-depsc','-r300','fig8.eps');
print(gcf,'-dpng','-r300','fig8.png');


%%

%% begin plotting
% figure('units','normalized','Position',[0.0 0.0 0.45, 1],'color','w');
nr=15;%number of stations in the first column
x0=0.1;y0=0.05;dy=0.13/2;dx=0;
%length: 0.5x0.5, 0.5x0.25, 0.25x0.5
%% axis XY
dh=(1-0.2)/3;dw=0.25;
dh1=0.04;%axis height
dy=0.04;

dh1=0.06;dy=dh1;
nr=3;
labels={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
ylm=[-999,1200];
ylm2=[-999,1200];
ylm3=[-999,1200];
% traces=[623:631];
traces=[50,80,250,480,520,630,750,900,910];

%trace:624,625,631: perfect
%311,

il=0;
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
for ir=1:3
%% first column
ix=traces(3*(ir-1)+1);
a1=axes('Parent',gcf,'Position',[x0,y0+dy+dh*(nr-ir),dw,dh1]);
plot(t,d11(:,ix),'k','linewidth',2); ylim(ylm);
if ir<=2
    set(gca,'xticklabel',[],'linewidth',2,'fontweight','bold');
else
    xlabel('Time (s)');set(gca,'linewidth',2,'fontweight','bold');
end
a1=axes('Parent',gcf,'Position',[x0,y0+2*dy+dh*(nr-ir),dw,dh1]);
plot(t,d1(:,ix),'k','linewidth',2); set(gca,'xticklabel',[],'linewidth',2,'fontweight','bold');ylim(ylm);
a1=axes('Parent',gcf,'Position',[x0,y0+3*dy+dh*(nr-ir),dw,dh1]);
plot(t,d2(:,ix),'k','linewidth',2); set(gca,'xticklabel',[],'linewidth',2,'fontweight','bold');ylim(ylm);
a1=axes('Parent',gcf,'Position',[x0,y0+4*dy+dh*(nr-ir),dw,dh1]);
plot(t,dn(:,ix),'k','linewidth',2); set(gca,'xticklabel',[],'linewidth',2,'fontweight','bold'); ylim(ylm);title(strcat('Channel:'," ",num2str(ix)));

%add component
a1=axes('Parent',gcf,'Position',[x0,y0+dy+dh*(nr-ir)+0.015,dw,dh1]);
text(0.01,dy*10,'SGK','color','r','Fontsize',10,'fontweight','bold');axis off;
a1=axes('Parent',gcf,'Position',[x0,y0+2*dy+dh*(nr-ir)+0.015,dw,dh1]);
text(0.01,dy*10,'KSVD','color','r','Fontsize',10,'fontweight','bold');axis off;
a1=axes('Parent',gcf,'Position',[x0,y0+3*dy+dh*(nr-ir)+0.015,dw,dh1]);
text(0.01,dy*10,'DCT','color','r','Fontsize',10,'fontweight','bold');axis off;
a1=axes('Parent',gcf,'Position',[x0,y0+4*dy+dh*(nr-ir)+0.015,dw,dh1]);
text(0.01,dy*10,'Raw data','color','r','Fontsize',10,'fontweight','bold');axis off;
%add label
il=il+1;
a1=axes('Parent',gcf,'Position',[x0,y0+5*dy+dh*(nr-ir)+0.015,dw,dh1]);
text(-0.1,0,labels(il),'color','k','Fontsize',15,'fontweight','bold');axis off;

%% Second column
ix=traces(3*(ir-1)+2);
a1=axes('Parent',gcf,'Position',[x0+0.3,y0+dy+dh*(nr-ir),dw,dh1]);
plot(t,d11(:,ix),'k','linewidth',2); ylim(ylm2);
if ir<=2
    set(gca,'xticklabel',[],'linewidth',2,'fontweight','bold');
else
    xlabel('Time (s)');set(gca,'linewidth',2,'fontweight','bold');
end
a1=axes('Parent',gcf,'Position',[x0+0.3,y0+2*dy+dh*(nr-ir),dw,dh1]);
plot(t,d1(:,ix),'k','linewidth',2); set(gca,'xticklabel',[],'linewidth',2,'fontweight','bold');ylim(ylm2);
a1=axes('Parent',gcf,'Position',[x0+0.3,y0+3*dy+dh*(nr-ir),dw,dh1]);
plot(t,d2(:,ix),'k','linewidth',2); set(gca,'xticklabel',[],'linewidth',2,'fontweight','bold');ylim(ylm2);
a1=axes('Parent',gcf,'Position',[x0+0.3,y0+4*dy+dh*(nr-ir),dw,dh1]);
plot(t,dn(:,ix),'k','linewidth',2); set(gca,'xticklabel',[],'linewidth',2,'fontweight','bold'); ylim(ylm2);title(strcat('Channel:'," ",num2str(ix)));

%add component
a1=axes('Parent',gcf,'Position',[x0+0.3,y0+dy+dh*(nr-ir)+0.015,dw,dh1]);
text(0.01,dy*10,'SGK','color','r','Fontsize',10,'fontweight','bold');axis off;
a1=axes('Parent',gcf,'Position',[x0+0.3,y0+2*dy+dh*(nr-ir)+0.015,dw,dh1]);
text(0.01,dy*10,'KSVD','color','r','Fontsize',10,'fontweight','bold');axis off;
a1=axes('Parent',gcf,'Position',[x0+0.3,y0+3*dy+dh*(nr-ir)+0.015,dw,dh1]);
text(0.01,dy*10,'DCT','color','r','Fontsize',10,'fontweight','bold');axis off;
a1=axes('Parent',gcf,'Position',[x0+0.3,y0+4*dy+dh*(nr-ir)+0.015,dw,dh1]);
text(0.01,dy*10,'Raw data','color','r','Fontsize',10,'fontweight','bold');axis off;
%label
il=il+1;
a1=axes('Parent',gcf,'Position',[x0+0.3,y0+5*dy+dh*(nr-ir)+0.015,dw,dh1]);
text(-0.1,0,labels(il),'color','k','Fontsize',15,'fontweight','bold');axis off;

%% Third column
ix=traces(3*(ir-1)+3);
a1=axes('Parent',gcf,'Position',[x0+0.6,y0+dy+dh*(nr-ir),dw,dh1]);
plot(t,d11(:,ix),'k','linewidth',2); ylim(ylm3);
if ir<=2
    set(gca,'xticklabel',[],'linewidth',2,'fontweight','bold');
else
    xlabel('Time (s)');set(gca,'linewidth',2,'fontweight','bold');
end
a1=axes('Parent',gcf,'Position',[x0+0.6,y0+2*dy+dh*(nr-ir),dw,dh1]);
plot(t,d1(:,ix),'k','linewidth',2); set(gca,'xticklabel',[],'linewidth',2,'fontweight','bold');ylim(ylm3);
a1=axes('Parent',gcf,'Position',[x0+0.6,y0+3*dy+dh*(nr-ir),dw,dh1]);
plot(t,d2(:,ix),'k','linewidth',2); set(gca,'xticklabel',[],'linewidth',2,'fontweight','bold');ylim(ylm3);
a1=axes('Parent',gcf,'Position',[x0+0.6,y0+4*dy+dh*(nr-ir),dw,dh1]);
plot(t,dn(:,ix),'k','linewidth',2); set(gca,'xticklabel',[],'linewidth',2,'fontweight','bold');ylim(ylm3);title(strcat('Channel:'," ",num2str(ix)));
%add component
a1=axes('Parent',gcf,'Position',[x0+0.6,y0+dy+dh*(nr-ir)+0.015,dw,dh1]);
text(0.01,dy*10,'SGK','color','r','Fontsize',10,'fontweight','bold');axis off;
a1=axes('Parent',gcf,'Position',[x0+0.6,y0+2*dy+dh*(nr-ir)+0.015,dw,dh1]);
text(0.01,dy*10,'KSVD','color','r','Fontsize',10,'fontweight','bold');axis off;
a1=axes('Parent',gcf,'Position',[x0+0.6,y0+3*dy+dh*(nr-ir)+0.015,dw,dh1]);
text(0.01,dy*10,'DCT','color','r','Fontsize',10,'fontweight','bold');axis off;
a1=axes('Parent',gcf,'Position',[x0+0.6,y0+4*dy+dh*(nr-ir)+0.015,dw,dh1]);
text(0.01,dy*10,'Raw data','color','r','Fontsize',10,'fontweight','bold');axis off;
%label
il=il+1;
a1=axes('Parent',gcf,'Position',[x0+0.6,y0+5*dy+dh*(nr-ir)+0.015,dw,dh1]);
text(-0.1,0,labels(il),'color','k','Fontsize',15,'fontweight','bold');axis off;

end
print(gcf,'-depsc','-r300','fig9.eps');
print(gcf,'-dpng','-r300','fig9.png');



%% fig10
figure('units','normalized','Position',[0.2 0.4 1 0.8]);
for ia=1:64
    subplot(8,8,ia);plot([1:l1],D(:,ia),'k','linewidth',2);
    
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
    axis off;
    xlim([1,l1]);
end
annotation(gcf,'rectangle',[0.126 0.101 0.782 0.834],'linewidth',2);
print(gcf,'-depsc','-r300','real_atom0_new.eps');

figure('units','normalized','Position',[0.2 0.4 0.5 1]);
for ia=1:64
    subplot(16,8,ia);plot([1:l1],Dksvd(:,ia),'k','linewidth',2);
    
    if ia==1
        ylim([-0.5,0.5]);
    end
        if ismember(ia,[12,13,14,18,53,54,57,58])
           subplot(16,8,ia);plot([1:l1],Dksvd(:,ia),'r','linewidth',2);
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
    axis off;
    xlim([1,l1]);
end
annotation(gcf,'rectangle',[0.126 0.101 0.782 0.834],'linewidth',2);
% print(gcf,'-depsc','-r300','real_atom1_new.eps');

% figure('units','normalized','Position',[0.2 0.4 1 0.8]);
for ia=1:64
    subplot(16,8,ia+64);plot([1:l1],Dsgk(:,ia),'k','linewidth',2);
    
    if ia==1
        ylim([-0.5,0.5]);
    end
        if ismember(ia,[6,7,8,9,10,11,12,13,14,32,47,52,54,55])
           subplot(16,8,ia+64);plot([1:l1],Dsgk(:,ia),'r','linewidth',2);
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
    axis off;
    xlim([1,l1]);
end
annotation(gcf,'rectangle',[0.126 0.101 0.782 0.834/2],'linewidth',2);
annotation(gcf,'rectangle',[0.126 0.101 0.782 0.834],'linewidth',2);

annotation(gcf,'textbox',...
    [0.07 0.8483 0.125 0.101],...
    'String','(a)','Fontsize',20,'fontweight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textbox',...
    [0.07 0.427 0.125 0.101],...
    'String','(b)','Fontsize',20,'fontweight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');
print(gcf,'-depsc','-r300','fig10.eps');
print(gcf,'-dpng','-r300','fig10.png');


