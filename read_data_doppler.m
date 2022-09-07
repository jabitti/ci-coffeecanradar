%MIT IAP Radar Course 2011
%Resource: Build a Small Radar System Capable of Sensing Range, Doppler, 
%and Synthetic Aperture Radar Imaging 
%
%Gregory L. Charvat

%Process Doppler vs. Time Intensity (DTI) plot

%NOTE: set Vtune to 3.2V to stay within ISM band and change fc to frequency
%below

clear all;
close all;

thresh = -45;
maxV = 40; % m/s
minV = 3; % m/s

%read the raw data .wave file here
[Y,FS] = audioread('perimeter_D_04.wav');
%constants
c = 3E8; %(m/s) speed of light

%radar parameters
Tp = 0.250; %(s) pulse time
N = Tp*FS; %# of samples per pulse
fc = 2295.27E6; %(Hz) Center frequency (connected VCO Vtune to +4.94)

% Create signal vector, s
s = Y(:,1);
clear Y;

%create doppler vs. time plot data set here
for ii = 1:round(size(s,1)/N)-1
    sif(ii,:) = s(1+(ii-1)*N:ii*N);
end

%subtract the average DC term here
sif = sif - mean(s);

zpad = 8*N/2;

%doppler vs. time plot:
v = dbv(ifft(sif,zpad,2));
v = v(:,1:size(v,2)/2);
mmax = max(max(v));
%calculate velocity
delta_f = linspace(0, FS/2, size(v,2)); %(Hz)
lambda=c/fc;
velocity = delta_f*lambda/2;
%calculate time
time = linspace(1,Tp*size(v,1),size(v,1)); %(sec)
%plot
vplot = v-mmax;
figure;
imagesc(velocity,time,vplot,[thresh, 0]);
colorbar;
xlim([3 40]); %limit velocity axis
xlabel('Velocity (m/sec)');
ylabel('time (sec)');

% create mesh plot
maxdelta = round(2*maxV/lambda);
meshdata = vplot(1:length(time),1:maxdelta);
figure;mesh(velocity(1:maxdelta),time, meshdata)


oneVec = ones([size(v,1),1]);
vpos = v-min(min(v));
meanplot = mean(vpos);
comp = oneVec * meanplot;
% figure;plot(velocity,meanplot)
vplot2 = v - comp;
vplot2 = vplot2 - max(max(vplot2));
% figure;imagesc(velocity,time,comp,[-35, 0])
% colorbar;
% xlim([0 40]); %limit velocity axis
% xlabel('Velocity (m/sec)');
% ylabel('time (sec)');figure;imagesc(velocity,time,v-comp,[-35, 0]);
figure;imagesc(velocity,time,vplot2,[thresh, 0])
colorbar;
xlim([0 40]); %limit velocity axis
xlabel('Velocity (m/sec)');
ylabel('time (sec)');
% figure;
% for ii = 1:size(sif,1)
%     subplot(2,1,1);
%     plot(sif(ii,:))
%     subplot(2,1,2);
%     plot(abs(fft(sif(ii,:))))
%     pause(0.3);
% end

% Doppler Movie
% vMax = 450;
% vplot3 = vplot2(:,1:vMax);
% vel3 = velocity(1:vMax);
% for ii = 1:size(vplot,1)
%     plot(vel3, vplot3(ii,:))
%     axis([0 vel3(vMax), -70 0])
%     pause(0.3);
% end

% Rainfall Plot
% figure;
% temp = -70 * ones(size(vplot2));
% for ii = 1:size(vplot2,1)
%     temp(1:ii,:) = vplot2(1:ii,:);
%     imagesc(velocity,time,temp,[-35, 0]);
%     colorbar;
%     xlim([0 40]); %limit velocity axis
%     xlabel('Velocity (m/sec)');
%     ylabel('time (sec)');
%     pause(Tp);
% end
