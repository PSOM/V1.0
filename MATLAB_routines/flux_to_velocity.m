function [u1,v1,w1] = flux_to_velocity(filename,xf,yf,zf,xc,yc)
% Function to convert from flux on the cell faces to velocity on the cell
% faces in PSOM. The flux on the cell face is divided by the area of the
% face. The top grid face moves due to variable sea surface height. This is
% accounted for in the calculation.
% MAF October 11, 2019

% Inputs
% filename: filename (as string) of 'face' file (e.g. 'face_022000.cdf'),
% which should have the fluxes on cell faces (uf, vf, wf) and the sea
% surface height (h)
% xf, yf, zf: coordinates of faces in x, y, and z directions

% Outputs
% Velocity on cell faces (u1, v1, w1)

h = ncread(filename,'h');
[X,Y] = meshgrid(xc,yc);

[~,dX,dZ] = meshgrid(1:321,diff(xf),diff(zf));
[X1,Y1] = meshgrid(xc(2:257),yf);
h1 = interp2(X,Y,h',X1,Y1);
dZ(:,:,end) = h1'-zf(end-1);
yarea = dX.*dZ;
vf = ncread(filename,'vf');
v1 = vf./yarea*1e5;

[dY,~,dZ] = meshgrid(diff(yf),1:257,diff(zf));
[X1,Y1] = meshgrid(xf,yc(2:321));
h1 = interp2(X,Y,h',X1,Y1);
dZ(:,:,end) = h1'-zf(end-1);
xarea = dY.*dZ;
uf = ncread(filename,'uf');
u1 = uf./xarea*1e5;

[dY,dX,~] = meshgrid(diff(yf),diff(xf),1:65);
zarea = dY.*dX;
wf = ncread(filename,'wf');
w1 = wf./zarea*1e4;