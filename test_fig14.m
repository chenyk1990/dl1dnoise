
clc;clear;close all;
load data/real3.mat

d=dd;

figure;
subplot(3,1,1);plot(dd(:,1));%E
subplot(3,1,2);plot(dd(:,2));%N
subplot(3,1,3);plot(dd(:,3));%Z
d=d(1:10:end,:);
d(:,1)=d(:,1:end-2);
d(:,2)=d(:,3:end);
d(:,3)=d(:,1:end-2);

figure;
subplot(3,1,1);plot(d(:,1));%E
subplot(3,1,2);plot(d(:,2));%N
subplot(3,1,3);plot(d(:,3));%Z

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
figure;
for ia=1:16
    subplot(4,4,ia);plot(dct(:,ia));
end

%% decompose the image into patches:
X=yc_patch(d,1,l1,1,l1/2,1);
fprintf('There are %d patches\n',size(X,2));

%% OMP using DCT
nd=size(X,2);
K=3;
ph=0.57;
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


%% K-SVD
param.T=K;      %sparsity level
param.D=D;    %initial D
param.niter=30; %number of K-SVD iterations to perform; default: 10
param.mode=1;   %1: sparsity; 0: error
param.K=c2;     %number of atoms, dictionary size
tic
[Dksvd,Gksvd]=yc_ksvd(X,param); 
toc

Gksvd0=Gksvd;
Gksvd=yc_pthresh(Gksvd0,'ph',ph);
X1=Dksvd*Gksvd;
[n1,n2]=size(d);
d1=yc_patch_inv(X1,1,n1,n2,l1,1,l1/2,1);

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

%% option2 (to make the signal leakage smallest, using an non-stationary strategy)
X2=yc_patch(d(775:end,:),1,l1,1,l1/2,1);
param.T=30;      %sparsity level
[Dsgk2,Gsgk2]=yc_sgk(X2,param); 
Gsgk2=yc_pthresh(Gsgk2,'ph',100); 
X22=Dsgk2*Gsgk2;
d11_2=yc_patch_inv(X22,1,1200-774,n2,l1,1,l1/2,1);
d11=[d11(1:774,:);d11_2];

%% plot
t=[0:1200-1]*0.1;
arr=768;
x1=60;x2=90;
figure('units','normalized','Position',[0.2 0.4 0.55, 0.9],'color','w');
subplot(6,3,1);plot(t,d(:,1),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Raw (E)','Fontsize',12,'fontweight','normal');ylabel('Amplitude (nm/s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');plot([x1,x1],[-10000,10000],'b','Linewidth',2);plot([x2,x2],[-10000,10000],'b','Linewidth',2);plot([x1,x2],[10000,10000],'b','Linewidth',2);plot([x1,x2],[-10000,-10000],'b','Linewidth',2);
subplot(6,3,4);plot(t,d(:,2),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Raw (N)','Fontsize',12,'fontweight','normal');ylabel('Amplitude (nm/s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');plot([x1,x1],[-10000,10000],'b','Linewidth',2);plot([x2,x2],[-10000,10000],'b','Linewidth',2);plot([x1,x2],[10000,10000],'b','Linewidth',2);plot([x1,x2],[-10000,-10000],'b','Linewidth',2);
subplot(6,3,7);plot(t,d(:,3),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Raw (Z)','Fontsize',12,'fontweight','normal');ylabel('Amplitude (nm/s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');plot([x1,x1],[-10000,10000],'b','Linewidth',2);plot([x2,x2],[-10000,10000],'b','Linewidth',2);plot([x1,x2],[10000,10000],'b','Linewidth',2);plot([x1,x2],[-10000,-10000],'b','Linewidth',2);
subplot(6,3,2);plot(t,d11(:,1),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Denoised (E)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');plot([x1,x1],[-10000,10000],'b','Linewidth',2);plot([x2,x2],[-10000,10000],'b','Linewidth',2);plot([x1,x2],[10000,10000],'b','Linewidth',2);plot([x1,x2],[-10000,-10000],'b','Linewidth',2);plot([x1,x1],[-10000,10000],'b','Linewidth',2);plot([x2,x2],[-10000,10000],'b','Linewidth',2);plot([x1,x2],[10000,10000],'b','Linewidth',2);plot([x1,x2],[-10000,-10000],'b','Linewidth',2);
subplot(6,3,5);plot(t,d11(:,2),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Denoised (N)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');plot([x1,x1],[-10000,10000],'b','Linewidth',2);plot([x2,x2],[-10000,10000],'b','Linewidth',2);plot([x1,x2],[10000,10000],'b','Linewidth',2);plot([x1,x2],[-10000,-10000],'b','Linewidth',2);
subplot(6,3,8);plot(t,d11(:,3),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Denoised (Z)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');plot([x1,x1],[-10000,10000],'b','Linewidth',2);plot([x2,x2],[-10000,10000],'b','Linewidth',2);plot([x1,x2],[10000,10000],'b','Linewidth',2);plot([x1,x2],[-10000,-10000],'b','Linewidth',2);
subplot(6,3,3);plot(t,d(:,1)-d11(:,1),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Noise (E)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');plot([x1,x1],[-10000,10000],'b','Linewidth',2);plot([x2,x2],[-10000,10000],'b','Linewidth',2);plot([x1,x2],[10000,10000],'b','Linewidth',2);plot([x1,x2],[-10000,-10000],'b','Linewidth',2);
subplot(6,3,6);plot(t,d(:,2)-d11(:,2),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Noise (N)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');plot([x1,x1],[-10000,10000],'b','Linewidth',2);plot([x2,x2],[-10000,10000],'b','Linewidth',2);plot([x1,x2],[10000,10000],'b','Linewidth',2);plot([x1,x2],[-10000,-10000],'b','Linewidth',2);
subplot(6,3,9);plot(t,d(:,3)-d11(:,3),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Noise (Z)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');plot([x1,x1],[-10000,10000],'b','Linewidth',2);plot([x2,x2],[-10000,10000],'b','Linewidth',2);plot([x1,x2],[10000,10000],'b','Linewidth',2);plot([x1,x2],[-10000,-10000],'b','Linewidth',2);

subplot(6,3,10);plot(t,d(:,1),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Zoom Raw (E)','Fontsize',12,'fontweight','normal');ylabel('Amplitude (nm/s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');xlim([x1,x2]);
subplot(6,3,13);plot(t,d(:,2),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Zoom Raw (N)','Fontsize',12,'fontweight','normal');ylabel('Amplitude (nm/s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');xlim([x1,x2]);
subplot(6,3,16);plot(t,d(:,3),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);xlabel('Time (s)','Fontsize',12,'fontweight','normal');title('Zoom Raw (Z)','Fontsize',12,'fontweight','normal');ylabel('Amplitude (nm/s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');xlim([x1,x2]);
subplot(6,3,11);plot(t,d11(:,1),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Zoom Denoised (E)','Fontsize',12,'fontweight','normal');ylabel('Amplitude (nm/s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');xlim([x1,x2]);xlim([x1,x2]);
subplot(6,3,14);plot(t,d11(:,2),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Zoom Denoised (N)','Fontsize',12,'fontweight','normal');ylabel('Amplitude (nm/s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');xlim([x1,x2]);
subplot(6,3,17);plot(t,d11(:,3),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);xlabel('Time (s)','Fontsize',12,'fontweight','normal');title('Zoom Denoised (Z)','Fontsize',12,'fontweight','normal');ylabel('Amplitude (nm/s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');xlim([x1,x2]);
subplot(6,3,12);plot(t,d(:,1)-d11(:,1),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Zoom Noise (E)','Fontsize',12,'fontweight','normal');ylabel('Amplitude (nm/s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');xlim([x1,x2]);
subplot(6,3,15);plot(t,d(:,2)-d11(:,2),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);title('Zoom Noise (N)','Fontsize',12,'fontweight','normal');ylabel('Amplitude (nm/s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');xlim([x1,x2]);
subplot(6,3,18);plot(t,d(:,3)-d11(:,3),'k','Linewidth',1.5);hold on;plot([arr,arr]*0.1,[-10000,10000],'r--','Linewidth',2);ylim([-10000,10000]);xlabel('Time (s)','Fontsize',12,'fontweight','normal');title('Zoom Noise (Z)','Fontsize',12,'fontweight','normal');xlabel('Time (s)','Fontsize',12,'fontweight','normal');set(gca,'Linewidth',1.5,'Fontsize',12,'Fontweight','normal');xlim([x1,x2]);

print(gcf,'-depsc','-r300','fig14.eps');
print(gcf,'-dpng','-r300','fig14.png');


