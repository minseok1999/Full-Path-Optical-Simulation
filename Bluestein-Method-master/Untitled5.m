clear all;
clc;

pixel0=8;
dim=625;
mx0=625;
my0=625;
%%%%% L0(mm)=diameter=8micrometerxdim(1)(mm)%%
L0=(dim)*pixel0;
%%%%everything written in units of micrometer%
[xx,yy]=meshgrid(-(mx0-1)/2:(mx0-1)/2,-(my0-1)/2:(my0-1)/2);
x0=xx.*pixel0;
y0=yy.*pixel0;
%No aperture%
waist=0.9e3;
gg=sign(1+sign(x0)).*sign(1+sign(y0));

%{
    % Grid in cylindrical coordinates:
    Rad = sqrt(x0.^2+y0.^2);
    Angle = angle(x0+1i.*y0)+pi; % Matrix with all angles starting left-center


gg=((-1+i).*GenModesLG([1 0], waist, Rad, Angle)+(1+i).*GenModesLG([-1 0], waist, Rad, Angle)+(0).*GenModesLG([0 0], waist, Rad, Angle))./2;
%}


%%%%%%%%%
%Input plot
%Input defining
h=abs(gg).^2; h=h./max(h(:)); %%%
%unit in mm in plotting, not micrometer%%
fig(1)=figure('Name','Input Gaussian amplitude');
surf(xx,yy,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'XTick',[-600:200:600],'YTick',[-600:200:600],'FontSize',24) 
%{
surf(xx.*L0./mx0./1e3,yy.*L0./my0./1e3,h), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[0,0.25,0.5,0.75,1]), shading interp;
set(gca,'XTick',[-4:1:4],'YTick',[-4:1:4],'FontSize',24)   
%phase plot%
H=angle(gg);%
fig(2)=figure('Name','Modulated phase');
surf(xx.*L0./mx0./1e3,yy.*L0./my0./1e3,H), colormap hot,axis equal, axis tight, view([0, -270]), colorbar('Ticks',[-1.5,0,1.5]), shading interp;
set(gca,'XTick',[-4:1:4],'YTick',[-4:1:4],'FontSize',24) 
%}