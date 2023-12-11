
%%%%%%Light: Science and Applications, DOI:10.1038/s41377-020-00362-z;
%%%%%%unit: um

clear all;
clc;
close all;
close all hidden;

matfile = 'unitary_trial';	
Hologram=load(fullfile('C:\Users\rmins\Desktop\MPLC\MPLC_system\WaveFrontMatching_Code\trials',matfile));	
Hologram=Hologram.Hologram;
Hologram=ones(1024,1024,5);
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

%%%%%%%%%%%%

waist=0.9e3;


    % Grid in cylindrical coordinates:
    Rad = sqrt(x0.^2+y0.^2);
    Angle = angle(x0+1i.*y0)+pi; % Matrix with all angles starting left-center

g=(GenModesLG([1 0], waist, Rad, Angle)+GenModesLG([-1 0], waist, Rad, Angle))./sqrt(2);

%%%%%%%%%

norm1=(sum(abs(g),'all')^2)




g=g.*exp(1i.*phi1);
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

tic;
d=800e3;                                 
% L=lamda*d/pixel0;                       
L1=L0;
x1start=-L1./2;
x1end=L1./2;
y1start=-L1./2;
y1end=L1./2; 
mx1=1024;
my1=1024;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g1,pixel1]=Scalar_Bluestein(g,mx0,my0,pixel0,d,x1start,x1end,y1start,y1end,mx1,my1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
x1=linspace(x1start,x1end,mx1);                
y1=linspace(y1start,y1end,my1);                
[x1,y1]=meshgrid(x1,y1);   
%Quickcheck
h=abs(g1).^2; h=h./max(h(:)); %%%
%%%%%%
fig(3)=figure('Name','1st amplitude');
surf(x1./1e3,y1./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(g1);%
fig(4)=figure('Name','1st phase');
surf(x1./1e3,y1./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%%%%%
g1=g1.*exp(1i.*Hologram(:,:,2));
%%%%%
L2=L0;
x2start=-L2./2;
x2end=L2./2;
y2start=-L2./2;
y2end=L2./2; 
mx2=1024;
my2=1024;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g2,pixel2]=Scalar_Bluestein(g1,mx1,my1,pixel1,d,x2start,x2end,y2start,y2end,mx2,my2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x2=linspace(x2start,x2end,mx2);                
y2=linspace(y2start,y2end,my2);                
[x2,y2]=meshgrid(x2,y2);     
%Quickcheck3
h=abs(g2).^2; h=h./max(h(:)); %%%
%%%%%%
fig(5)=figure('Name','2nd amplitude');
surf(x2./1e3,y2./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(g2);%
fig(6)=figure('Name','2nd phase');
surf(x2./1e3,y2./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%%%%%
g2=g2.*exp(1i.*Hologram(:,:,3));

%check
h=abs(g2).^2; h=h./max(h(:)); %%%
%%%%%%
fig(7)=figure('Name','final amplitude');
surf(x2./1e3,y2./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(g2);%
fig(8)=figure('Name','final phase');
surf(x2./1e3,y2./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%%%%%

norm2=(sum(abs(g2),'all')^2)

%%%%%%%

L3=L1;
x3start=-L3./2;
x3end=L3./2;
y3start=-L3./2;
y3end=L3./2; 
mx3=1024;
my3=1024;
d=d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g3,pixel3]=Scalar_Bluestein(g2,mx2,my2,pixel2,d,x3start,x3end,y3start,y3end,mx3,my3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x3=linspace(x3start,x3end,mx3);                
y3=linspace(y3start,y3end,my3);                
[x3,y3]=meshgrid(x3,y3);    

%%%mode overlap with desired beam

 % Grid in cylindrical coordinates:
    Rad = sqrt(x3.^2+y3.^2);
    Angle = angle(x3+1i.*y3)+pi; % Matrix with all angles starting left-center
    
gdot=((i).*GenModesLG([1 0], waist, Rad, Angle)+(1+2i).*GenModesLG([-1 0], waist, Rad, Angle)+(1+i).*GenModesLG([0 0], waist, Rad, Angle))./sqrt(8);
overlap=(abs(sum(g3.*conj(gdot),'all'))^2)./(sum(abs(g3).^2,'all'))./(sum(abs(gdot).^2,'all'))

%check
h=abs(g3).^2; h=h./max(h(:)); %%%
%%%%%%
fig(9)=figure('Name','far field amplitude');
surf(x3./1e3,y3./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(g3);%
fig(10)=figure('Name','far field phase');
surf(x3./1e3,y3./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%%%%%

norm3=(sum(abs(g3),'all')^2)

h=abs(gdot).^2; h=h./max(h(:)); %%%
%%%%%%
fig(11)=figure('Name','desired amplitude');
surf(x3./1e3,y3./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(g3);%
fig(12)=figure('Name','desired phase');
surf(x3./1e3,y3./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%%%%%


norm4=(sum(abs(gdot),'all')^2)


