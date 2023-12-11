%%%%%%Scalar calculation
%%%%%%Light: Science and Applications, DOI:10.1038/s41377-020-00362-z;
%%%%%%unit: um

%%%%%automated
clear all;
clc;
close all;
close all hidden;

  
tic;
global lamda k
lamda=780e-3;
k=2*pi/lamda;
%sampling point interval in micrometer=pixel0%
pixel0=8;
%phase modulation%
%%%%%%
phi=fliplr(readmatrix('001L1qutritIntensity'));	
%%%%%%
toc

tic;
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
w0=1.5e3;  z0=w0^2*pi./lamda; R=f*(1+(z0/f)^2);%%
g=exp(-(x0.^2+y0.^2)./w0^2).*exp(1i*phi);
%Input plot
%Input defining
h=abs(g).^2; h=h./max(h(:)); %%%
%unit in mm in plotting, not micrometer%%
fig(1)=figure('Name','Input Gaussian amplitude');
surf(xx.*L0./mx0./1e3,yy.*L0./my0./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(g);%
fig(2)=figure('Name','Modulated phase');
surf(xx.*L0./mx0./1e3,yy.*L0./my0./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%Input plot done%%%
%%%%projection angle theta
theta=0.35/180*pi;
diff=exp(1i.*k.*sin(theta).*(x0+L0/2));

toc
%%%%%%
g=g.*diff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

d=f;                                 
% L=lamda*d/pixel0;                       
L1=L0;
x1start=-L1./2;
x1end=L1./2;
y1start=-L1./2;
y1end=L1./2; 
mx1=1920;
my1=1920;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g1,pixel1]=Scalar_Bluestein(g,mx0,my0,pixel0,d,x1start,x1end,y1start,y1end,mx1,my1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
x1=linspace(x1start,x1end,mx1);                
y1=linspace(y1start,y1end,my1);                
[x1,y1]=meshgrid(x1,y1);         
%%%%%%%
lens=g1.*exp(-1i.*k./2./f.*(x1.^2+y1.^2));
%%%to the focus%%%
L2=L0./20;
x2start=-L2./2;
x2end=L2./2;
y2start=-L2./2;
y2end=L2./2; 
mx2=1920;
my2=1920;
d2=f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g2,pixel2]=Scalar_Bluestein(lens,mx1,my1,pixel1,d2,x2start,x2end,y2start,y2end,mx2,my2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x2=linspace(x2start,x2end,mx2);                
y2=linspace(y2start,y2end,my2);                
[x2,y2]=meshgrid(x2,y2);     
%%%Cutting with pinhole%%%
pinhole=sign(1-sign(x2.^2+y2.^2-(L2./5).^2));
A0=pinhole;
g2=A0.*g2;
%Quickcheck3
h=abs(g2).^2; h=h./max(h(:)); %%%
%%%%%%
fig(3)=figure('Name','1st pinhole amplitude');
surf(x2./1e3,y2./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:0.1:4],'XTick',[-4:0.1:4],'FontSize',24)   
%phase plot%
H=angle(g2);%
fig(4)=figure('Name','1st pinhole phase');
surf(x2./1e3,y2./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:0.1:4],'XTick',[-4:0.1:4],'FontSize',24) 
%%%%%%%
%%%free propagation to lens%%%
L3=L1;
x3start=-L3./2;
x3end=L3./2;
y3start=-L3./2;
y3end=L3./2; 
mx3=1920;
my3=1920;
d3=f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g3,pixel3]=Scalar_Bluestein(g2,mx2,my2,pixel2,d3,x3start,x3end,y3start,y3end,mx3,my3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x3=linspace(x3start,x3end,mx3);                
y3=linspace(y3start,y3end,my3);                
[x3,y3]=meshgrid(x3,y3);    
x3=x3.*L3./L1;
y3=y3.*L3./L1;
g3=g3.*exp(-1i.*k./2./f.*(x3.^2+y3.^2));
%%%%%
%%%propagation after the 2nd lens%%%
L4=L1;
x4start=-L4./2;
x4end=L4./2;
y4start=-L4./2;
y4end=L4./2; 
mx4=1920;
my4=1920;
d4=f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g4,pixel4]=Scalar_Bluestein(g3,mx3,my3,pixel3,d3,x4start,x4end,y4start,y4end,mx4,my4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x4=linspace(x4start,x4end,mx4);                
y4=linspace(y4start,y4end,my4);                
[x4,y4]=meshgrid(x4,y4);     
toc

%Quickcheck
h=abs(g4).^2; h=h./max(h(:)); %%%
%%%%%%
fig(5)=figure('Name','at 400 mm amplitude');
surf(x4./1e3,y4./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(g4);%
fig(6)=figure('Name','at 400mm phase');
surf(x4./1e3,y4./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%%%%%

%%%%2nd SLM%%%
tic;
%phase modulation%
phi2=fliplr(readmatrix('L1qutrit_projection'));
dim2=size(phi2);
mx4=dim2(1);
my4=dim2(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g4=g4.*exp(1i*phi2);


f=125e3;
%%%%%
d5=f;                                                   
L5=L1;
x5start=-L5./2;
x5end=L5./2;
y5start=-L5./2;
y5end=L5./2; 
mx5=1920;
my5=1920;
%%%%propagation to the lens%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g5,pixel5]=Scalar_Bluestein(g4,mx4,my4,pixel4,d5,x5start,x5end,y5start,y5end,mx5,my5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%Lensing%%%%
x5=linspace(x5start,x5end,mx5);                
y5=linspace(y5start,y5end,my5);                
[x5,y5]=meshgrid(x5,y5);    
x5=x5.*L5./L1;
y5=y5.*L5./L1;
lens=g5.*exp(-1i.*k./2./f.*(x5.^2+y5.^2));
%%%to the focus%%%
L6=L0./40;
x6start=-L6./2;
x6end=L6./2;
y6start=-L6./2;
y6end=L6./2; 
mx6=1920;
my6=1920;
d6=f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g6,pixel6]=Scalar_Bluestein(lens,mx5,my5,pixel5,d6,x6start,x6end,y6start,y6end,mx6,my6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x6=linspace(x6start,x6end,mx6);                
y6=linspace(y6start,y6end,my6);                
[x6,y6]=meshgrid(x6,y6);     
%%%Cutting with pinhole%%%
pinhole=sign(1-sign(x6.^2+y6.^2-(L6./5).^2));
A0=pinhole;
g6=A0.*g6;
%%%%
h=abs(g6).^2; h=h./max(h(:)); %%%
fig(7)=figure('Name','after cutting amplitude');
surf(x6./1e3,y6./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:0.1:4],'XTick',[-4:0.1:4],'FontSize',24)   
%phase plot%
H=angle(g6);%
fig(8)=figure('Name','after cutting phase');
surf(x6./1e3,y6./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:0.1:4],'XTick',[-4:0.1:4],'FontSize',24) 
%%%
%%%free propagation to lens%%%
L7=L1;
x7start=-L7./2;
x7end=L7./2;
y7start=-L7./2;
y7end=L7./2; 
mx7=1920;
my7=1920;
d7=f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g7,pixel7]=Scalar_Bluestein(g6,mx6,my6,pixel6,d7,x7start,x7end,y7start,y7end,mx7,my7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x7=linspace(x7start,x7end,mx7);                
y7=linspace(y7start,y7end,my7);                
[x7,y7]=meshgrid(x7,y7);    
x7=x7.*L7./L1;
y7=y7.*L7./L1;
g7=g7.*exp(-1i.*k./2./f.*(x7.^2+y7.^2));
%%%%%
%%%propagation after the 4th lens%%%
L8=L1;
x8start=-L8./2;
x8end=L8./2;
y8start=-L8./2;
y8end=L8./2; 
mx8=1920;
my8=1920;
d8=f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g8,pixel8]=Scalar_Bluestein(g7,mx7,my7,pixel7,d7,x8start,x8end,y8start,y8end,mx8,my8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x8=linspace(x8start,x8end,mx8);                
y8=linspace(y8start,y8end,my8);                
[x8,y8]=meshgrid(x8,y8);     
toc

%%%%%Tight focusing

f=8e3;
lens=g8.*exp(-1i.*k./2./f.*(x8.^2+y8.^2));
%%%to the focus%%%
L9=L0./500;
x9start=-L9./2;
x9end=L9./2;
y9start=-L9./2;
y9end=L9./2; 
mx9=1920;
my9=1920;
d9=f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g9,pixel9]=Scalar_Bluestein(lens,mx8,my8,pixel8,d9,x9start,x9end,y9start,y9end,mx9,my9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x9=linspace(x9start,x9end,mx9);                
y9=linspace(y9start,y9end,my9);                
[x9,y9]=meshgrid(x9,y9);     


%%%%
h=abs(g9).^2; h=h./max(h(:)); %%%
fig(9)=figure('Name','at coupler amplitude');
surf(x9./1e3,y9./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:0.01:4],'XTick',[-4:0.01:4],'FontSize',24)   
%phase plot%
H=angle(g9);%
fig(10)=figure('Name','at coupler phase');
surf(x9./1e3,y9./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:0.01:4],'XTick',[-4:0.01:4],'FontSize',24) 

%%%%mode overlap%%%%%%
%%%%diameter of single mode fiber is 5micrometer

gdot=exp(-(x9.^2+y9.^2)./(5)^2);
overlap=(abs(sum(g9.*conj(gdot),'all'))^2)./(sum(abs(g9).^2,'all'))./(sum(abs(gdot).^2,'all'))
