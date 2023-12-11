%%%%%%Scalar calculation
%%%%%%Light: Science and Applications, DOI:10.1038/s41377-020-00362-z;
%%%%%%unit: um
%%%%%%Making CSV

clear all;
clc;
tic;

global lamda k
lamda=780e-3;
k=2*pi/lamda;
%sampling point interval in micrometer=pixel0%
pixel0=4;
%phase modulation%
phi=readmatrix('0012L12qutrit.csv');
dim=size(phi);
mx0=dim(1);
my0=dim(2);
%%%%% L0(mm)=diameter=4micrometerxdim(1)(mm)%%
L0=(dim(1))*pixel0;
%%%%everything written in units of micrometer%
[xx,yy]=meshgrid(-(mx0-1)/2:(mx0-1)/2,-(my0-1)/2:(my0-1)/2);
x0=xx.*pixel0;
y0=yy.*pixel0;
%No aperture%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=100e3;%
%define w0, waist function at z=0
w0=1.66e3;  z0=w0^2*pi./lamda; R=f*(1+(z0/f)^2);%%
g=exp(-(x0.^2+y0.^2)./w0^2).*exp(1i*phi);
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
toc
%Input plot
%Input defining
h=abs(g).^2; h=h./max(h(:)); %%%

%%%Input plot done%%%
x1=linspace(x1start,x1end,mx1);                
y1=linspace(y1start,y1end,my1);                
[x1,y1]=meshgrid(x1,y1);                   
lens=g1.*exp(-1i.*k./2./f.*(x1.^2+y1.^2));
%%%to the focus%%%
tic;
L2=L1;
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
x2=linspace(x2start,x2end,mx2);                
y2=linspace(y2start,y2end,my2);                
[x2,y2]=meshgrid(x2,y2);     
toc
%%%Cutting with pinhole%%%
tic;
pinhole=sign(1-sign(x2.^2+y2.^2-(L2./10).^2));
A0=pinhole;
g2=A0.*g2;
%%%free propagation to lens%%%
L3=L1;
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
x3=x3.*L3./L1;
y3=y3.*L3./L1;
g3=g3.*exp(-1i.*k./2./f.*(x3.^2+y3.^2));
toc
%%%%%
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

%%%%2nd SLM%%%

%phase modulation%
phi2=readmatrix('L12qutrit_projection.csv');
dim2=size(phi2);
mx4=dim2(1);
my4=dim2(2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g4=g4.*exp(1i*phi2);

%%%%%
d5=f;                                                   
L5=L1;
x5start=-L5./2;
x5end=L5./2;
y5start=-L5./2;
y5end=L5./2; 
mx5=1080;
my5=1080;
%%%%propagation to the lens%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g5,pixel5]=Scalar_Bluestein(g4,mx4,my4,pixel4,d5,x5start,x5end,y5start,y5end,mx5,my5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
toc
%%%%Lensing%%%%
x5=linspace(x5start,x5end,mx5);                
y5=linspace(y5start,y5end,my5);                
[x5,y5]=meshgrid(x5,y5);    
x5=x5.*L5./L1;
y5=y5.*L5./L1;
lens=g5.*exp(-1i.*k./2./f.*(x5.^2+y5.^2));
%%%to the focus%%%
tic;
L6=L0./20;
x6start=-L6./2;
x6end=L6./2;
y6start=-L6./2;
y6end=L6./2; 
mx6=1080;
my6=1080;
d6=f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g6,pixel6]=Scalar_Bluestein(lens,mx5,my5,pixel5,d6,x6start,x6end,y6start,y6end,mx6,my6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x6=linspace(x6start,x6end,mx6);                
y6=linspace(y6start,y6end,my6);                
[x6,y6]=meshgrid(x6,y6);     
toc

%%%%
writematrix(g6,fullfile('C:\Users\rmins\Desktop\array','output12.csv'));
