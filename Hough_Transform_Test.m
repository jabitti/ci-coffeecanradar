% Only need to go from 0-180 degrees since they repeat after that
close all; clear; clc;
points = diag(ones(1,10));

transform = hough(points);
imagesc(transform);
mesh(transform);