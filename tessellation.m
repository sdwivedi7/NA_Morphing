function [p,t] = tessellation(file_name)
data = importdata(file_name);
%%
x = data(:,1);
y = data(:,2);
minX = min(abs(x - [x(2:end);x(1)]));
minY = min(abs(y - [y(2:end);y(1)]));
p = 10*max(minX,minY);
pts = [-1.2:p:1.2];
[xmesh,ymesh] = meshgrid(pts);
in = inpolygon(xmesh,ymesh,x,y);
%%
xmesh = xmesh(in); ymesh = ymesh(in);
x = [x;xmesh];
y = [y;ymesh];
%%
t = delaunay(x,y);
%%
p = [x';y'];
t = t';
% x = data(:,1);
% y = data(:,2);
end