% Only need to go from 0-180 degrees since they repeat after that
close all; clear; clc;
points = diag(ones(1,10));
%points = diag(ones(1,20));
%points(end,:) = [];
points = rot90(points);

transform = hough(points);

[maxRows, iRows] = max(transform);
[maxColumns, deg] = max(maxRows);
deg = deg - 90;
[point1,point2] = size(points);

[rhoMax, degMax] = size(transform);
rho = point1 * (iRows(deg)/(rhoMax-1));
% iColumn gives the rho
%rho = 5*sqrt(2);

x_0 = rho * cosd(deg);
y_0 = rho * sind(deg);
m1 = -(tand(-deg))^-1;
x = 1:1:10;

m2 = -(1/tand(-deg));
b = rho/sind(deg);

line1 = m1 .* (x - x_0) + y_0;
line2 = m2 .* x + b;

figure(); hold on;
%plot(1:10,points,'color','b');
plot(x,line1,'color','r'); axis([0,10 0,22]);
plot(x,line2,'color','b');
%imagesc(line2);

%imagesc(transform);
figure(); mesh(transform);
%figure(); imagesc(points);