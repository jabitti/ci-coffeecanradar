% Only need to go from 0-180 degrees since they repeat after that
close all; clear; clc;
points = diag(ones(1,10));
%points = rot90(points);

transform = hough(points);

[maxRows, iRows] = max(transform);
[maxColumns, deg] = max(maxRows);

[rhoMax, degMax] = size(transform);
rho = iRows(deg);
% iColumn gives the rho

x_0 = rho * cosd(deg);
y_0 = rho * sind(deg);
m1 = -(tand(deg))^-1;
x = 0:0.01:10;

line = m1 .* (x - x_0) + y_0;

figure(); plot(x,line); axis([0,10 0,22]);

%imagesc(transform);
figure(); mesh(transform);
figure(); imagesc(points);