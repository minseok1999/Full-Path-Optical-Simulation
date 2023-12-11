clear all;
clc;
close all;
close all hidden;

tic;
global lamda k
lamda=808e-3;
k=2*pi/lamda;
%sampling point interval in micrometer=pixel0%
pixel0=8;
%%%%%%
mx0=1080;
my0=1080;
%%%%% L0(mm)=diameter=8micrometerxdim(1)(mm)%%
L0=mx0*pixel0;
%%%%everything written in units of micrometer%
[xx,yy]=meshgrid(-(mx0-1)/2:(mx0-1)/2,-(my0-1)/2:(my0-1)/2);
x0=xx.*pixel0;
y0=yy.*pixel0;


ModeTypes=[1,0];
waist=0.9e3;

toc

    % Grid in cylindrical coordinates:
    Rad = sqrt(x0.^2+y0.^2);
    Angle = angle(x0+1i.*y0)+pi; % Matrix with all angles starting left-center

    
g=GenModesLG(ModeTypes, waist, Rad, Angle);

%Input plot
%Input defining
h=abs(g).^2; h=h./max(h(:)); %%%
%unit in mm in plotting, not micrometer%%
fig(1)=figure('Name','amplitude');
surf(x0./1e3,y0./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(g);%
fig(2)=figure('Name','phase');
surf(x0./1e3,y0./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'YTick',[-4:1:4],'XTick',[-4:1:4],'FontSize',24) 
%%%Input plot done%%%


