%%%%%%Scalar calculation of the converging spherical wave
%%%%%%Light: Science and Applications, DOI:10.1038/s41377-020-00362-z;
%%%%%%unit: um

clear all;
clc;
tic;

global lamda k
lamda=800e-3;
k=2*pi/lamda;

my0=1081;
mx0=1081;
pixel0=8;
L0=(mx0-1)*pixel0;
[xx,yy]=meshgrid(-(my0-1)/2:(my0-1)/2,-(my0-1)/2:(my0-1)/2);
x0=xx.*pixel0;
y0=yy.*pixel0;
Aperture=sign(1-sign(xx.^2+yy.^2-((my0-1)./2).^2));
A0=Aperture;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=600e3;%

%define w0, waist function at z=0
w0=1e3;  z0=w0^2*pi./lamda; R=f*(1+(z0/f)^2);%%
%Phase
phi=atan(y0./x0);%

g=A0.*exp(-1i.*k./2./f.*(x0.^2+y0.^2)).*exp(-1i.*k./2.*(x0.^2+y0.^2)./R).*exp(-(x0.^2+y0.^2)./w0^2);
% figure
% imshow(angle(g),[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=f;                                 
% L=lamda*d/pixel0;                       
L=0.2e3;
x1start=-L./2;
x1end=L./2;
y1start=-L./2;
y1end=L./2; 
mx1=1081;
my1=1081;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g1,pixel1]=Scalar_Bluestein(g,mx0,my0,pixel0,d,x1start,x1end,y1start,y1end,mx1,my1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=abs(g1).^2;
P=angle(g1);
x1=linspace(x1start,x1end,mx1);                
y1=linspace(y1start,y1end,my1);                
[x1,y1]=meshgrid(x1,y1); 
toc
%Input plot
%Input defining
h=abs(g).^2; h=h./max(h(:)); %%%
%
figure('Name','Input amplitude')
surf(xx./1e3,yy./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-0.1:0.1:0.1],'XTick',[-0.1:0.1:0.1],'FontSize',24)   
%Input phase%
H=angle(g);%
%phase plot%
figure('Name','Input phase')
surf(xx./1e3,yy./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-0.1:0.1:0.1],'XTick',[-0.1:0.1:0.1],'FontSize',24) 
%%Input plot done%%
figure('Name','Output amplitude')
surf(x1./1e3,y1./1e3,I), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-0.1:0.1:0.1],'XTick',[-0.1:0.1:0.1],'FontSize',24)   
% saveas(gcf,'Intensity.png');
figure('Name','Output phase')
surf(x1./1e3,y1./1e3,P), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-0.1:0.1:0.1],'XTick',[-0.1:0.1:0.1],'FontSize',24) 
