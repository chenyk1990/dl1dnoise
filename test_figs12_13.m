% Script to plot Figure 10
% BY Yangkang Chen
% Dec, 23, 2021

clc;clear;close all;


addpath(genpath('./subroutines/'));

%% load data
eq=zeros(2000,960);
[n1,n2]=size(eq);
ii=3;
% if ~ismember(ii,[14,16,17,27,47,52])
%     load(strcat('/Users/chenyk/dasdenoising/mat_raw/eq-',num2str(ii),'.mat'));
% end
load data/real2.mat
eq=d1(:,:);
eq=d1;
d_bp=yc_bandpass(d1,0.0005,0,200,6,6,0,0);%
% d_bpfk=d_bp-yc_fk_dip(d_bp,0.02);%
d_bpfk=d_bp;
d=d_bpfk;
figure;yc_imagesc([eq,d,eq-d]);

% load(strcat('/Users/chenyk/dasdenoising/mat_bpsomffk/eq-',num2str(ii),'.mat'));
% load data/real11.mat
% d=d1;
% figure;yc_imagesc([eq,d,eq-d]);

eq=eq(:,100:800);
d=d(:,100:800);
%% denoising
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
figure;
for ia=1:16
    subplot(4,4,ia);plot(dct(:,ia));
end

% decompose the image into patches:
X=yc_patch(d,1,l1,1,l1/2,1);


% OMP using DCT
nd=size(X,2);
K=3;
ph=2;
tic
for i2=1:nd
    G(:,i2)=yc_omp0(D,X(:,i2),K);
end
toc

%further constrain it to be sparser
G=yc_pthresh(G,'ph',ph);
X2=D*G;

[n1,n2]=size(d);
d2=yc_patch_inv(X2,1,n1,n2,l1,1,l1/2,1);
figure;yc_imagesc([d,d2,d-d2]);

% figure('units','normalized');
% imagesc(G);colormap(jet);colorbar;caxis([-0.5,0.5]);%colorbar;
% ylabel('Atom NO','Fontsize',16);
% xlabel('Patch NO','Fontsize',16);
% title('Coefficients Matrix','Fontsize',16);
% set(gca,'Linewidth',1.5,'Fontsize',16);


% %% K-SVD
% param.T=K;      %sparsity level
% param.D=D;    %initial D
% param.niter=30; %number of K-SVD iterations to perform; default: 10
% param.mode=1;   %1: sparsity; 0: error
% param.K=c2;     %number of atoms, dictionary size
% tic
% [Dksvd,Gksvd]=yc_ksvd(X,param); 
% toc
% 
% Gksvd0=Gksvd;
% Gksvd=yc_pthresh(Gksvd0,'ph',ph);
% X1=Dksvd*Gksvd;
% [n1,n2]=size(d);
% d1=yc_patch_inv(X1,1,n1,n2,l1,1,l1/2,1);

%% SGK
param.T=K;      %sparsity level
param.D=D;    %initial D
param.niter=30; %number of K-SVD iterations to perform; default: 10
param.mode=1;   %1: sparsity; 0: error
param.K=c2;     %number of atoms, dictionary size
tic
[Dsgk,Gsgk]=yc_sgk(X,param); 
toc
Gsgk0=Gsgk;
Gsgk=yc_pthresh(Gsgk0,'ph',ph);
X11=Dsgk*Gsgk;
d11=yc_patch_inv(X11,1,n1,n2,l1,1,l1/2,1);

% save fig10.mat 

%%cost
% DCT:15.16s
% SGK:35.00s

% NOs=[1,3,20,10,25,11,2];
% labels={...                                          %P-arrival sample NO from the SEGY file
%     'FORGE\_78-32\_iDASv3-P11\_UTC190423150554.sgy',... %24169
%     'FORGE\_78-32\_iDASv3-P11\_UTC190423213209.sgy',... 
%     'FORGE\_78-32\_iDASv3-P11\_UTC190426070723.sgy',... %24811
%     'FORGE\_78-32\_iDASv3-P11\_UTC190426062208.sgy',... %26090
%     'FORGE\_78-32\_iDASv3-P11\_UTC190426110008.sgy',... %4921
%     'FORGE\_78-32\_iDASv3-P11\_UTC190426062553.sgy',... %8934
%     'FORGE\_78-32\_iDASv3-P11\_UTC190423182409.sgy'};   %4210

% eq=zeros(2000,960);
[n1,n2]=size(d);
t=[0:n1]*0.0005;
ngap=50;
x=1:n2*5+4*ngap;
%% first one
ii=1;
indt1=10:300;indx1=300:700;

d1_z1=d2(indt1,indx1);
d1_z2=d11(indt1,indx1);
comp1=[eq,zeros(n1,ngap),d2,zeros(n1,ngap),eq-d2,zeros(n1,ngap),d11,zeros(n1,ngap),eq-d11]; 



comp2=comp1;
%% combined figure
figure('units','normalized','Position',[0.0 0.0 0.6, 0.5],'color','w');
yc_imagesc(comp1(:,:),98,1,x,t(:));
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(n2/2,-0.02,'Raw data','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text((n2*2+ngap*3)/2+n2,-0.02,'DCT','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text((n2*2+ngap*3)/2+n2*3+ngap*2,-0.02,'SGK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.02,'(a)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
% text(50,0.46,labels{1},'color','b','Fontsize',14,'fontweight','bold','HorizontalAlignment','left');
hold on;
plot([indx1(1)+n2+ngap,indx1(1)+n2+ngap],t([indt1(1),indt1(end)]),'r','linewidth',2);
plot([indx1(end)+n2+ngap,indx1(end)+n2+ngap],t([indt1(1),indt1(end)]),'r','linewidth',2);
plot([indx1(1)+n2+ngap,indx1(end)+n2+ngap],t([indt1(1),indt1(1)]),'r','linewidth',2);
plot([indx1(1)+n2+ngap,indx1(end)+n2+ngap],t([indt1(end),indt1(end)]),'r','linewidth',2);

plot([indx1(1)+n2*3+ngap*3,indx1(1)+n2*3+ngap*3],t([indt1(1),indt1(end)]),'r','linewidth',2);
plot([indx1(end)+n2*3+ngap*3,indx1(end)+n2*3+ngap*3],t([indt1(1),indt1(end)]),'r','linewidth',2);
plot([indx1(1)+n2*3+ngap*3,indx1(end)+n2*3+ngap*3],t([indt1(1),indt1(1)]),'r','linewidth',2);
plot([indx1(1)+n2*3+ngap*3,indx1(end)+n2*3+ngap*3],t([indt1(end),indt1(end)]),'r','linewidth',2);

annotation(gcf,'arrow',[0.394 0.364],[0.806 0.674],'linewidth',2,'color','r');
annotation(gcf,'textarrow',[0.377 0.377],[0.621 0.591],'String','Artifacts','linewidth',2,'color','g','fontsize',15,'fontweight','bold');
annotation(gcf,'arrow',[0.71 0.68],[0.806 0.674],'linewidth',2,'color','r');
annotation(gcf,'arrow',[0.697 0.697],[0.621 0.591],'linewidth',2,'color','g');
annotation(gcf,'textarrow',[0.5023 0.4849],[0.7056 0.79333],'String',{'Signal leakage'},'linewidth',2,'fontsize',15,'fontweight','bold','color','g');

% subplot(2,1,2);yc_imagesc(comp2(1:1000,:),98,1,x,t(1:1000));
% ylabel('Time (s)','Fontsize',14,'fontweight','bold');
% xlabel('Channel','Fontsize',14,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
% text(n2/2,-0.02,'Raw data','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text((n2*2+ngap*3)/2+n2,-0.02,'BP+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text((n2*2+ngap*3)/2+n2*3+ngap*2,-0.02,'BP+SOMF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.02,'(b)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
% text(50,0.46,labels{2},'color','b','Fontsize',14,'fontweight','bold','HorizontalAlignment','left');
% hold on;
% plot([indx2(1)+n2+ngap,indx2(1)+n2+ngap],t([indt2(1),indt2(end)]),'r','linewidth',2);
% plot([indx2(end)+n2+ngap,indx2(end)+n2+ngap],t([indt2(1),indt2(end)]),'r','linewidth',2);
% plot([indx2(1)+n2+ngap,indx2(end)+n2+ngap],t([indt2(1),indt2(1)]),'r','linewidth',2);
% plot([indx2(1)+n2+ngap,indx2(end)+n2+ngap],t([indt2(end),indt2(end)]),'r','linewidth',2);
% 
% plot([indx2(1)+n2*3+ngap*3,indx2(1)+n2*3+ngap*3],t([indt2(1),indt2(end)]),'r','linewidth',2);
% plot([indx2(end)+n2*3+ngap*3,indx2(end)+n2*3+ngap*3],t([indt2(1),indt2(end)]),'r','linewidth',2);
% plot([indx2(1)+n2*3+ngap*3,indx2(end)+n2*3+ngap*3],t([indt2(1),indt2(1)]),'r','linewidth',2);
% plot([indx2(1)+n2*3+ngap*3,indx2(end)+n2*3+ngap*3],t([indt2(end),indt2(end)]),'r','linewidth',2);
% 
% annotation(gcf,'arrow',[0.394 0.364],[0.371 0.332],'linewidth',2,'color','r');
% annotation(gcf,'textarrow',[0.307 0.307],[0.281 0.311],'String','Artifacts','linewidth',2,'color','g','fontsize',15,'fontweight','bold');
% annotation(gcf,'arrow',[0.71 0.68],[0.371 0.332],'linewidth',2,'color','r');
% annotation(gcf,'arrow',[0.627 0.627],[0.281 0.311],'linewidth',2,'color','g');


%% add zooming framebox

a1=axes('Parent',gcf,'Position',[0.287,0.38,0.149,0.3]);
yc_imagesc(d1_z1,70,2);axis off;

a1=axes('Parent',gcf,'Position',[0.600,0.38,0.149,0.3]);
yc_imagesc(d1_z2,70,2);axis off;


% a1=axes('Parent',gcf,'Position',[0.287,0.18,0.149,0.15]);
% yc_imagesc(d2_z1,100,2);axis off;
% 
% a1=axes('Parent',gcf,'Position',[0.600,0.18,0.149,0.15]);
% yc_imagesc(d2_z2,100,2);axis off;

print(gcf,'-depsc','-r300','fig12.eps');



inds=20:20:n2;
traces=d11(:,inds);
traces0=eq(:,inds);

dn=traces;
d0=traces0;
nsta=30;nlta=80;
[ O,R ] = dl_picker_stalta(dn,nsta, nlta);
[ O0,R0 ] = dl_picker_stalta(d0,nsta, nlta);

times0=[O0-1]*0.0005;
times=[O-1]*0.0005;
%name
for ii=1:30
    stname{ii}=strcat('Channel:'," ",num2str(inds(ii)));
end

%% begin plotting
figure('units','normalized','Position',[0.0 0.0 0.45, 1],'color','w');
nr=7;%number of stations in the first column
x0=0.1;y0=0.1;dy=0;dx=0;
%length: 0.5x0.5, 0.5x0.25, 0.25x0.5
%% axis XY
dh=(1-0.2)/nr;dw=0.37;
dh1=dh;%axis height
t=[0:1999]*0.0005;
for ir=nr:-1:1
    wav0=traces0(:,ir);
    wav=traces(:,ir);
    a1=axes('Parent',gcf,'Position',[x0,y0+dy+dh*(nr-ir),dw,dh1]);
    plot(t,wav0,'k','linewidth',2); hold on; axis off;
    plot(t,wav,'r','linewidth',2);
    plot([times0(ir),times0(ir)],[min(wav),max(wav)],'g','linewidth',2);
    plot([times(ir),times(ir)],[min(wav),max(wav)],'b','linewidth',2);
    
    wav0=traces0(:,ir+15);
    wav=traces(:,ir+15);
    a1=axes('Parent',gcf,'Position',[x0+0.5,y0+dy+dh*(nr-ir),dw,dh1]);
    plot(t,wav0,'k','linewidth',2);hold on; axis off; 
    plot(t,wav,'r','linewidth',2);
    plot([times0(ir+15),times0(ir+15)],[min(wav),max(wav)],'g','linewidth',2);
    plot([times(ir+15),times(ir+15)],[min(wav),max(wav)],'b','linewidth',2);
    
end
legend('Raw waveform','Denoised waveform','Picked arrival from raw data','Picked arrival from denoised data','Position',[x0+0.15,y0-0.1,0.6,0.1],'NumColumns',4);
legend('boxoff');
% 
%% add station name
for ir=nr:-1:1
a1=axes('Parent',gcf,'Position',[0.02,y0+dh*(nr-ir)+dh/2+0.01,dw,dh1]);
text(-0.035,0,stname{ir},'color','k','Fontsize',10,'fontweight','bold');axis off;

a1=axes('Parent',gcf,'Position',[0.02+0.5,y0+dh*(nr-ir)+dh/2+0.01,dw,dh1]);
text(-0.035,0,stname{ir+15},'color','k','Fontsize',10,'fontweight','bold');axis off;
end
%
%
%% add source info
dw2=(1-x0)/5.0;
a1=axes('Parent',gcf,'Position',[0,0.93,1,dh1]);
text(0.5,0,'Microseismic event detection','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');axis off;
print(gcf,'-depsc','-r300','fig13.eps');



