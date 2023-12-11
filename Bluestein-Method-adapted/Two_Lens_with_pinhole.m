%%%%%%Scalar calculation
%%%%%%Light: Science and Applications, DOI:10.1038/s41377-020-00362-z;
%%%%%%unit: um

clear all;
clc;
tic;

global lamda k
lamda=780e-3;
k=2*pi/lamda;
%sampling point interval in micrometer=pixel0%
pixel0=4;
%phase modulation%
phi=readmatrix('004L4qutrit.csv');
dim=size(phi);
mx0=dim(1);
my0=dim(2);
%%%%% L0(mm)=diameter=8micrometerxdim(1)(mm)%%
L0=(dim(1))*pixel0;
%%%%everything written in units of micrometer%
[xx,yy]=meshgrid(-(mx0-1)/2:(mx0-1)/2,-(my0-1)/2:(my0-1)/2);
x0=xx.*pixel0;
y0=yy.*pixel0;
%No aperture%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=100e3;%
%define w0, waist function at z=0
w0=0.75e3;  z0=w0^2*pi./lamda; R=f*(1+(z0/f)^2);%%
g=exp(-(x0.^2+y0.^2)./w0^2).*exp(-(x0.^2+y0.^2)./w0^2).*exp(1i*phi);
% figure
% imshow(angle(g),[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=f;                                 
% L=lamda*d/pixel0;                       
L1=L0;
x1start=-L1./2;
x1end=L1./2;
y1start=-L1./2;
y1end=L1./2; 
mx1=1080;
my1=1080;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g1,pixel1]=Scalar_Bluestein(g,mx0,my0,pixel0,d,x1start,x1end,y1start,y1end,mx1,my1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
n1=pixel1
toc
%Input plot
%Input defining
h=abs(g).^2; h=h./max(h(:)); %%%
%unit in mm in plotting, not micrometer%%
figure('Name','Input Gaussian amplitude')
surf(xx.*L0./mx0./1e3,yy.*L0./my0./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(g);%
figure('Name','Modulated phase')
surf(xx.*L0./mx0./1e3,yy.*L0./my0./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%Input plot done%%%
x1=linspace(x1start,x1end,mx1);                
y1=linspace(y1start,y1end,my1);                
[x1,y1]=meshgrid(x1,y1);    
lens=g1.*exp(-1i.*k./2./f.*(x1.^2+y1.^2));
I=abs(lens).^2;I=I./max(I(:));
P=angle(lens);
%after the lens%
figure('Name','At Lens amplitude')
surf(x1./L1.*mx1.*L1./mx1./1e3,y1./1e3,I), colormap hot,axis equal,axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
% saveas(gcf,'Intensity.png');
figure('Name','At Lens phase')
surf(x1./1e3,y1./1e3,P), colormap hot,axis equal,axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%%%%
%%%to the focus%%%
tic;
L2=L1./10;
x2start=-L2./2;
x2end=L2./2;
y2start=-L2./2;
y2end=L2./2; 
mx2=1080;
my2=1080;
d2=f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g2,pixel2]=Scalar_Bluestein(lens,mx1,my1,pixel1,d2,x2start,x2end,y2start,y2end,mx2,my2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n2=pixel2
x2=linspace(x2start,x2end,mx2);                
y2=linspace(y2start,y2end,my2);                
[x2,y2]=meshgrid(x2,y2);     
toc
%%%%
h=abs(g2).^2; h=h./max(h(:)); %%%
%
figure('Name','before pinhole amplitude')
surf(x2.*L2./mx2./1e3,y2.*L2./mx2./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:0.02:4],'XTick',[-4:0.02:4],'FontSize',24)   
%phase plot%
H=angle(g2);%
figure('Name','before pinhole phase')
surf(x2.*L2./mx2./1e3,y2.*L2./mx2./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:0.02:4],'XTick',[-4:0.02:4],'FontSize',24) 
%%%
%%%Cutting with pinhole%%%
tic;
pinhole=sign(1-sign(x2.^2+y2.^2-(L2./5).^2));
A0=pinhole;
g2=A0.*g2;
%%%Plot again after cutting %%%%%
%%%%
h=abs(g2).^2; h=h./max(h(:)); %%%
%
figure('Name','after pinhole amplitude')
surf(x2.*L2./mx2./1e3,y2.*L2./mx2./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:0.02:4],'XTick',[-4:0.02:4],'FontSize',24)   
%phase plot%
H=angle(g2);%
figure('Name','after pinhole phase')
surf(x2.*L2./mx2./1e3,y2.*L2./mx2./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:0.02:4],'XTick',[-4:0.02:4],'FontSize',24) 
%%%%%
%%%free propagation to lens%%%
L3=L1./3;
x3start=-L3./2;
x3end=L3./2;
y3start=-L3./2;
y3end=L3./2; 
mx3=1080;
my3=1080;
d3=f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g3,pixel3]=Scalar_Bluestein(g2,mx2,my2,pixel2,d3,x3start,x3end,y3start,y3end,mx3,my3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x3=linspace(x3start,x3end,mx3);                
y3=linspace(y3start,y3end,my3);                
[x3,y3]=meshgrid(x3,y3);    
g3=g3.*exp(-1i.*k./2./f.*(x3.^2+y3.^2));
toc
%%%%%
h=abs(g3).^2; h=h./max(h(:)); %%%
figure('Name','2nd lens amplitude')
surf(x3.*L3./mx3./1e3,y3.*L3./my3./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:0.2:4],'XTick',[-4:0.2:4],'FontSize',24)   
%phase plot%
H=angle(g3);%
figure('Name','2nd lens phase')
surf(x3.*L3./mx3./1e3,y3.*L3./my3./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:0.2:4],'XTick',[-4:0.2:4],'FontSize',24) 

%%%propagation after the 2nd lens%%%
tic;
L4=L1;
x4start=-L4./2;
x4end=L4./2;
y4start=-L4./2;
y4end=L4./2; 
mx4=1080;
my4=1080;
d4=f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g4,pixel4]=Scalar_Bluestein(g3,mx3,my3,pixel3,d3,x4start,x4end,y4start,y4end,mx4,my4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x4=linspace(x4start,x4end,mx4);                
y4=linspace(y4start,y4end,my4);                
[x4,y4]=meshgrid(x4,y4);     
toc
%%%%
h=abs(g4).^2; h=h./max(h(:)); %%%
%
figure('Name','final plane amplitude')
surf(x4.*L4./mx4./1e3,y4.*L4./mx4./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(g4);%
figure('Name','final plane phase')
surf(x4.*L4./mx4./1e3,y4.*L4./mx4./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%