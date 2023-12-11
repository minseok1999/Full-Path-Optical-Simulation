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
%Phase%
mod=pi.*sign(1-sign(x0));
phi=atan(y0./x0)+mod;%
g=A0.*(x0.^2+y0.^2).^(1/2).*exp(-1i.*k./2./f.*(x0.^2+y0.^2)).*exp(-1i.*k./2.*(x0.^2+y0.^2)./R).*exp(-(x0.^2+y0.^2)./w0^2).*exp(-1i*phi);
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
figure('Name','Input amplitude')
f(1)=surf(xx.*8.64./1080,yy.*8.64./1080,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%Input phase%
H=angle(g);%
%phase plot%
figure('Name','Input phase')
f(2)=surf(xx.*8.64./1080,yy.*8.64./1080,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%Input plot done%%
%output spacing is 8.64e3/1080 /10e3 %
%output%
figure('Name','Output amplitude')
f(3)=surf(x1./1e3,y1./1e3,I), colormap hot,axis equal,axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
% saveas(gcf,'Intensity.png');
figure('Name','Output phase')
f(4)=surf(x1./1e3,y1./1e3,P), colormap hot,axis equal,axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 

for k=1:4
  
  saveas(f(k), fullfile('C:\Users\rmins\Desktop\matfig',sprintf('fig_%d',k)), 'jpeg')
 
end

