close all; clear; clc;

% Generate Array of Points
points = diag(ones(1,1000));
%points = diag(ones(1,20));
%points(end,:) = [];
%points = rot90(points);
[rows,cols] = size(points); % Get Size of Points Array For Later Scaling

% Take the Hough Transform, Save in "transform"
transform = hough(points);

% Extract the Highest Peak Rho and Theta
[maxRows, iRows] = max(transform);
[maxColumns, deg] = max(maxRows);
[rhoMax, ~] = size(transform); % Get Size of Transform for Scaling
rho = iRows(deg) - (rhoMax-1)/2.0; % iRows(deg) is Row of Maximum Value
    % [(rhoMax-1)/2] is the Scaled Maximum
deg = deg - 90; % Convert [1 to 180] to [-89 to 90]

m = -1/tand(deg);
r = 1:10;
c = m*r + rho/sind(deg);

for i = 1:rows
    j = m*i + rho/sind(deg);
    if j < 1
        j = 1;
    end
    points(i,round(j)) = 2;
end

%% Graphs
figure(); mesh(points);
figure(); plot(c,r); xlabel('c'); ylabel('r');
figure(); imagesc(points);