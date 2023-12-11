%%%%%%Scalar calculation of the converging spherical wave
%%%%%%Light: Science and Applications, DOI:10.1038/s41377-020-00362-z;
%%%%%%unit: um

clear all;
clc;
tic;

global lamda k
lamda=800e-3;
k=2*pi/lamda;
%1081-1=1080 for evenness, in an attempt to avoid singularity at sampling point(0,0)*
my0=1080;
mx0=1080;
%sampling point interval in micrometer=pixel0%
pixel0=8;
%%%%% L0(mm)=diameter=8micrometerx1080=8.64mm%%
L0=(mx0-1)*pixel0*1e3;
%%%%1e3 multiplied, everything written in units of micrometer%
[xx,yy]=meshgrid(-(my0-1)/2:(my0-1)/2,-(my0-1)/2:(my0-1)/2);
x0=xx.*pixel0;
y0=yy.*pixel0;
%circular aperture%
Aperture=sign(1-sign(xx.^2+yy.^2-((my0-1)./2).^2));
A0=Aperture;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=600e3;%
%define w0, waist function at z=0
w0=1e3;  z0=w0^2*pi./lamda; R=f*(1+(z0/f)^2);%%
%phase modulation%
phi=pi.*sign(1-sign(x0));
g=exp(-(x0.^2+y0.^2)./w0^2).*exp(-(x0.^2+y0.^2)./w0^2).*exp(1i*phi);
% figure
% imshow(angle(g),[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=f;                                 
% L=lamda*d/pixel0;                       
L=8.64e3;
x1start=-L./2;
x1end=L./2;
y1start=-L./2;
y1end=L./2; 
mx1=1080;
my1=1080;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g1,pixel1]=Scalar_Bluestein(g,mx0,my0,pixel0,d,x1start,x1end,y1start,y1end,mx1,my1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=abs(g1).^2;I=I./max(I(:));
P=angle(g1);
x1=linspace(x1start,x1end,mx1);                
y1=linspace(y1start,y1end,my1);                
[x1,y1]=meshgrid(x1,y1);     
toc
%Input plot
%Input defining
%Input spacing is set to be equal to output, i.e, by multiplying and dividing to yield 1*8.64/1080%
h=abs(g).^2; h=h./max(h(:)); %%%
%
figure('Name','Input Gaussian amplitude')
surf(xx.*8.64./1080,yy.*8.64./1080,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%Input phase%
H=angle(g);%
%phase plot%
figure('Name','Modulated phase')
surf(xx.*8.64./1080,yy.*8.64./1080,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%Input plot done%%
%output spacing is 8.64e3/1080 /10e3 %
%output%
figure('Name','At Lens amplitude')
surf(x1./1e3,y1./1e3,I), colormap hot,axis equal,axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
% saveas(gcf,'Intensity.png');
figure('Name','At Lens phase')
surf(x1./1e3,y1./1e3,P), colormap hot,axis equal,axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 

%after the lens%
tic;
lens=g1.*exp(-1i.*k./2./f.*(x1.^2+y1.^2));
L2=8.64e3;
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
aftl=abs(g2).^2;
aftl=aftl./max(aftl(:));
aftP=angle(g2);
x2=linspace(x2start,x2end,mx2);                
y2=linspace(y2start,y2end,my2);                
[x2,y2]=meshgrid(x2,y2);     
toc
%%%%
h=abs(g2).^2; h=h./max(h(:)); %%%
%
figure('Name','After amplitude')
surf(xx.*8.64./mx2,yy.*8.64./my2,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%Input phase%
H=angle(g2);%
%phase plot%
figure('Name','After phase')
surf(xx.*8.64./mx2,yy.*8.64./my2,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
phasemap;
