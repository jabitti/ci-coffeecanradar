%MIT IAP Radar Course 2011
%Resource: Build a Small Radar System Capable of Sensing Range, Doppler, 
%and Synthetic Aperture Radar Imaging 
%
%Gregory L. Charvat

%Process Range vs. Time Intensity (RTI) plot

%NOTE: set up-ramp sweep from 2-3.2V to stay within ISM band
%change fstart and fstop bellow when in ISM band

clear all;
close all;
clc;

rangeMovie = 'no';
dBthresh = -55;
minRg = 20;

collectionsfiles = dir('./Collections/*.wav');

%for i = 1:height(collectionsfiles)
for i = 3
    close all;
    %read the raw data .wav file here
    % Y contains an N x 2 matrix with the following:
    %   Column 1: radar return IF signal
    %   Column 2: sync square wave from function generator
    % FS contains the sample frequency of the recording software
    %   Typically 44100 samples/sec
    strcat("./Collections/",collectionsfiles(i).name);
    [Y,FS] = audioread(strcat("./Collections/",collectionsfiles(i).name));

    %constants
    c = 3E8; %(m/s) speed of light

    %radar parameters
    %%%%% ===== Time Between Vertical Steps
    Tp = 20E-3; %(s) pulse time
    N = Tp*FS; %# of samples per pulse
    fstart = 2225.8E6; %(Hz) LFM start frequency from performance pred. spreadsheet
    fstop = 2364.7E6; %(Hz) LFM stop frequency from performance pred. spreadsheet

    BW = fstop-fstart; %(Hz) transmtit bandwidth = difference of start and stop freqs
    f = linspace(fstart, fstop, N/2); %instantaneous transmit frequency

    %range resolution
    rr = c/(2*BW);
    %fprintf("\n%0.3f rr\n",rr);
    %%%%% ===== MAX RANGE
    % max_range
    max_range = rr*N/2;
    
    % Create separate vectors for the IF signal and the sync signal
    trig = Y(:,2);
    s = Y(:,1);
    clear Y;

    %parse the data here by triggering off rising edge of sync pulse
    count = 0;
    thresh = 0;
    start = (trig > thresh);
    for ii = 100:(size(start,1)-N)
        if start(ii) == 1 & mean(start(ii-11:ii-1)) == 0
            count = count + 1;
            sif(count,:) = s(ii:ii+N-1);
            time(count) = ii*1/FS;
        end
    end
    %check to see if triggering works
    % Ni = 1;
    % Nf = 15000;
    % plot(trig(Ni:Nf),'.b');
    % hold on;
    % plot(start(Ni:Nf),'.r');
    % hold off;
    % grid on;

    %subtract the average
    ave = mean(sif,1);
    for ii = 1:size(sif,1)
        sif(ii,:) = sif(ii,:) - ave;
    end

    zpad = 8*N/2;
    R = linspace(0,max_range,zpad/2);

    %RTI plot
    fig10 = figure(10);
    fig10.Renderer = 'Painters'; % this ensures that your figure will be vector/ps
    v = dbv(ifft(sif,zpad,2));
    S = v(:,1:size(v,2)/2);
    m = max(max(v));
    S1 = S-m; % RTI only
    imagesc(R,time,S1,[dBthresh, 0]);
    colorbar;
    ylabel('time (s)');
    xlabel('range (m)');
    title('RTI without clutter rejection');
    figpath = strcat("./CollectionResults/"+collectionsfiles(i).name(1:2)+ ...
        "-OG"+collectionsfiles(i).name(3:end-4));
    %%%print(figpath,'-dpng','-r600') % or change this to 600 (dpi) for crispier figs

    %2 pulse cancelor RTI plot
    fig20 = figure(20);
    fig20.Renderer = 'Painters'; % this ensures that your figure will be vector/ps
    sif2 = sif(2:size(sif,1),:)-sif(1:size(sif,1)-1,:);
    v = ifft(sif2,zpad,2);
    S=v;
    % for ii = 1:size(S,1)
    %     S(ii,:) = S(ii,:).*R.^(3/2); %Optional: magnitude scale to range
    % end
    S = dbv(S(:,1:size(v,2)/2));
    m = max(max(S));
    S2 = S-m; % RTI with 2 pulse cancellation
    imagesc(R,time,S2,[dBthresh, 0]);
    colorbar;
    ylabel('time (s)');
    xlabel('range (m)');
    title('RTI with 2-pulse cancelor clutter rejection');
    figpath = strcat("./CollectionResults/"+collectionsfiles(i).name(1:2)+...
        "-PC"+collectionsfiles(i).name(3:end-4));
    %%%print(figpath,'-dpng','-r600') % or change this to 600 (dpi) for crispier figs
    
%% Hough Transform
    %%% Isolate Main Target
    % S is the 3D graph of the return signal
    %figure(); mesh(S); title("OG Signal");
    % Plot the transpose of S (gives max of each row)
    %figure(); plot(max(S')); title("Max of Rows");
    
%     Smn = mean(S2);
%     S2norm = S2-Smn;
%     threshold = 20;
%     figure;mesh(S2norm);
%     S2bin = (S2norm > threshold) .* ones(size(S2norm));
%     figure(); imagesc(S2bin);
    
    [rows,cols] = size(S2);
    transform = hough(S2);
    [maxRows, iRows] = max(transform);
    [maxColumns, deg] = max(maxRows);
    [rhoMax, ~] = size(transform);
    rho = iRows(deg) - (rhoMax-1)/2.0;
        transform(iRows(deg+90),deg+90) = -50;
    deg = deg - 90;
    m = -1/tand(deg);
    for i = 1:rows
        j = m*i + rho/sind(deg);
        if j < 1
            j = 1;
        end
        S2(i,round(j)) = -10;
    end
    
    % Attempt to Extract Another?
%     [maxRows, iRows] = max(transform);
%     [maxColumns, deg] = max(maxRows);
%     rho = iRows(deg) - (rhoMax-1)/2.0;
%         transform(iRows(deg+90),deg+90) = -50;
%     deg = deg - 90;
%     m = -1/tand(deg);
%     for i = 1:rows
%         j = m*i + rho/sind(deg);
%         if j < 1
%             j = 1;
%         end
%         S(i,round(j)) = -10;
%     end
    figure(); mesh(transform);
    figure(); imagesc(S2);
    
    
%% Velocity Calculation

    velocity = ( (m*rows)-(m*0 + rho/sind(deg)) ) / ( rows*Tp );
    fprintf('Velocity: %0.4f m/s\n',velocity);
    fprintf('Velocity: %0.4f mph\n\n',velocity/2.23694);

    % figure(25);
    % plotminR = 20; % meters
    % plotmaxR = 120; % meters
    % delR = max_range/zpad;
    % Rzoom = linspace(plotminR,plotmaxR,floor((plotmaxR-plotminR)/delR));
    % imagesc(Rzoom,time,S2(:,floor(plotminR/delR):floor(plotmaxR/delR)),[dBthresh, 0]);
    % colorbar;
    % ylabel('time (s)');
    % xlabel('range (m)');
    % title('Zoomed RTI with 2-pulse cancelor clutter rejection');
    i = i + 1;
end

%% Plot single time slice
% figure;
% [Rw,Cl]=size(S0);
% width = 150;
% figure
% for ii=1:Rw
%     plot(S0(ii,1:width))
%     axis([0,width,-150,0])
%     pause(0.04)
% end
% %2 pulse mag only cancelor
% figure(30);
% clear v;
% for ii = 1:size(sif,1)-1
%     v1 = abs(ifft(sif(ii,:),zpad));
%     v2 = abs(ifft(sif(ii+1,:),zpad));
%     v(ii,:) = v2-v1;
% end
% S=v;
% R = linspace(0,max_range,zpad);
% for ii = 1:size(S,1)
%     S(ii,:) = S(ii,:).*R.^(3/2); %Optional: magnitude scale to range
% end
% S = dbv(S(:,1:size(v,2)/2));
% m = max(max(S));
% imagesc(R,time,S-m,[-20, 0]);
% colorbar;
% ylabel('time (s)');
% xlabel('range (m)');
% title('RTI with 2-pulse mag only cancelor clutter rejection');