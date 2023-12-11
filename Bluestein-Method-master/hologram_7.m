
%%%%%%Light: Science and Applications, DOI:10.1038/s41377-020-00362-z;
%%%%%%unit: um

clear all;
clc;
close all;
close all hidden;

Hologram=load('7holos00to10.mat');
Hologram=Hologram.Hologram;
tic;
global lamda k
lamda=808e-3;
k=2*pi/lamda;
%sampling point interval in micrometer=pixel0%
pixel0=8;
%phase modulation%
%%%%%%	
phi1=Hologram(:,:,1);
%%%%%%
toc

dim=size(phi1);
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
%define w0, waist function at z=0
w0=0.9e3;  
g=exp(-(x0.^2+y0.^2)./w0^2).*exp(1i.*phi1);
%Input plot
%Input defining
h=abs(g).^2; h=h./max(h(:)); %%%
%unit in mm in plotting, not micrometer%%
fig(1)=figure('Name','Input Gaussian amplitude');
surf(xx.*L0./mx0./1e3,yy.*L0./my0./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(g);%
fig(2)=figure('Name','Input phase');
surf(xx.*L0./mx0./1e3,yy.*L0./my0./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%Input plot done%%%
%%%Phase modulation
g=g.*exp(1i.*phi1);
%%%%%
tic;
d=800e3;                                 
% L=lamda*d/pixel0;                       
L1=L0;
x1start=-L1./2;
x1end=L1./2;
y1start=-L1./2;
y1end=L1./2; 
mx1=625;
my1=625;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g1,pixel1]=Scalar_Bluestein(g,mx0,my0,pixel0,d,x1start,x1end,y1start,y1end,mx1,my1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
x1=linspace(x1start,x1end,mx1);                
y1=linspace(y1start,y1end,my1);                
[x1,y1]=meshgrid(x1,y1);   

%%%%%%%
g1=g1.*exp(1i.*Hologram(:,:,2));
%%%%%
L2=L0;
x2start=-L2./2;
x2end=L2./2;
y2start=-L2./2;
y2end=L2./2; 
mx2=625;
my2=625;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g2,pixel2]=Scalar_Bluestein(g1,mx1,my1,pixel1,d,x2start,x2end,y2start,y2end,mx2,my2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x2=linspace(x2start,x2end,mx2);                
y2=linspace(y2start,y2end,my2);                
[x2,y2]=meshgrid(x2,y2);     

%%%%%%%
g2=g2.*exp(1i.*Hologram(:,:,3));

%%%%%%%
L3=L1;
x3start=-L3./2;
x3end=L3./2;
y3start=-L3./2;
y3end=L3./2; 
mx3=625;
my3=625;
d=d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g3,pixel3]=Scalar_Bluestein(g2,mx2,my2,pixel2,d,x3start,x3end,y3start,y3end,mx3,my3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x3=linspace(x3start,x3end,mx3);                
y3=linspace(y3start,y3end,my3);                
[x3,y3]=meshgrid(x3,y3);    

%%%%
g3=g3.*exp(1i.*Hologram(:,:,4));
%%%%%%%
L4=L1;
x4start=-L4./2;
x4end=L4./2;
y4start=-L4./2;
y4end=L4./2; 
mx4=625;
my4=625;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g4,pixel4]=Scalar_Bluestein(g3,mx3,my3,pixel3,d,x4start,x4end,y4start,y4end,mx4,my4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x4=linspace(x4start,x4end,mx4);                
y4=linspace(y4start,y4end,my4);                
[x4,y4]=meshgrid(x4,y4);    


%%%%
g4=g4.*exp(1i.*Hologram(:,:,5));
%%%%%%%
L5=L1;
x5start=-L5./2;
x5end=L5./2;
y5start=-L5./2;
y5end=L5./2; 
mx5=625;
my5=625;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g5,pixel5]=Scalar_Bluestein(g4,mx4,my4,pixel4,d,x5start,x5end,y5start,y5end,mx5,my5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x5=linspace(x5start,x5end,mx5);                
y5=linspace(y5start,y5end,my5);                
[x5,y5]=meshgrid(x5,y5);    

%%%%
g5=g5.*exp(1i.*Hologram(:,:,6));
%%%%%%%
L6=L1;
x6start=-L6./2;
x6end=L6./2;
y6start=-L6./2;
y6end=L6./2; 
mx6=625;
my6=625;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g6,pixel6]=Scalar_Bluestein(g5,mx5,my5,pixel5,d,x6start,x6end,y6start,y6end,mx6,my6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x6=linspace(x6start,x6end,mx6);                
y6=linspace(y6start,y6end,my6);                
[x6,y6]=meshgrid(x6,y6);    

%%%%
g6=g6.*exp(1i.*Hologram(:,:,7));

%check
h=abs(g6).^2; h=h./max(h(:)); %%%
%%%%%%
fig(3)=figure('Name','7th hologram amplitude');
surf(x6./1e3,y6./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(g6);%
fig(4)=figure('Name','7th hologram phase');
surf(x6./1e3,y6./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%%%%%

%%%%far field

%%%%%%%
L7=5.*L1;
x7start=-L7./2;
x7end=L7./2;
y7start=-L7./2;
y7end=L7./2; 
mx7=625;
my7=625;
d=10*d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g7,pixel7]=Scalar_Bluestein(g6,mx6,my6,pixel6,d,x7start,x7end,y7start,y7end,mx7,my7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x7=linspace(x7start,x7end,mx7);                
y7=linspace(y7start,y7end,my7);                
[x7,y7]=meshgrid(x7,y7);   


%check
h=abs(g7).^2; h=h./max(h(:)); %%%
%%%%%%
fig(5)=figure('Name','far field amplitude');
surf(x7./1e3,y7./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-40:5:40],'XTick',[-40:5:40],'FontSize',24)   
%phase plot%
H=angle(g7);%
fig(6)=figure('Name','far field phase');
surf(x7./1e3,y7./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-40:5:40],'XTick',[-40:5:40],'FontSize',24) 
%%%%%%%



for k=1:6
  
 saveas(fig(k), fullfile('C:\Users\rmins\Desktop\hologram',sprintf('00to10_%d',k)), 'jpeg')
  
end
