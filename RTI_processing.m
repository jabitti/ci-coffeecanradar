close all
rangeMovie = 'yes';
plotS2r = 'no';
plotS2L = 'yes';
plotS2mean = 'yes';
plotS2meanpoly = 'no';
plotS2rS2L = 'no';
trackTgt = 'no';
thresh = -50;
CFARthresh = 15;
%% RTI Plot Data
% S2: RTI 2-pulse clutter rejection data

%% Plot mean removed imagesc
S2r = S2-mean(S2);
S2r = S2r-max(max(S2r));
if strcmp(plotS2r,'yes')
    figure(10);
    imagesc(R,time,S2r,[thresh, 0]);
    colorbar;
    ylabel('time (s)');
    xlabel('range (m)');
    title('RTI with mean range profile removed');
end

%% PLot smoothed mean removed
S2mean = mean(S2);
% smooth S2mean
lpf_max=100;
S2meanF = fft(S2mean);
S2meanLPF = zeros([1,length(S2meanF)]);
S2meanLPF(1:lpf_max) = S2meanF(1:lpf_max);
S2meanSm = fliplr(abs(ifft(S2meanLPF)));

% least squares fit of S2meansm to 3rd order polynomial
x = 1:length(S2meanSm);
X = [x'.^6 x'.^5 x'.^4 x'.^3 x'.^2 x' ones([length(x) 1])];
coeff = inv(X'*X)*X'*S2meanSm';
figure(20);plot(R,S2meanSm,R,(X*coeff)','r-')
xlabel('Range (R) [m]');
ylabel('Amplitude (s) [dB]');
title('RTI with smoothed mean range profile removed')

S2L = S2-S2meanSm;
S2L = S2L-max(max(S2L));
if strcmp(plotS2L, 'yes')
    figure(21);
    plot(R,abs(S2meanSm))
    xlabel('Range (R) [m]');
    ylabel('Amplitude (s) [dB]');
    title('RTI with smoothed mean range profile removed + Polyfit')
    figure(22);
    imagesc(R,time,S2L,[thresh,0]);
    colorbar;
    title('RTI with smoothed mean range profile removed')
    ylabel('time (s)');
    xlabel('range (m)');
    
    figure(24);
    mesh(R(1:150),time(1:end-1),S2L(:,1:150));
    title('RTI with smoothed mean range profile removed')
    ylabel('time (s)');
    xlabel('range (m)');
    colorbar;
end

if strcmp(plotS2rS2L, 'yes')
    figure(30)
    subplot(1,2,1)
    imagesc(R,time,S2r,[thresh, 0]);
    colorbar;
    ylabel('time (s)');
    xlabel('range (m)');
    title('RTI with mean range profile removed');
    
    subplot(1,2,2)
    imagesc(R,time,S2L,[thresh, 0]);
    colorbar;
    ylabel('time (s)');
    xlabel('range (m)');
    title('RTI with smoothed mean range profile removed');
end
%% CFAR Processing
[rows,cols]=size(S2L);
nrg = mean(S2L,2);
nrg2D = nrg*ones([1 cols]);
S2LC = S2L-nrg2D;
S2LC = (S2LC>CFARthresh).*S2LC;
figure(30);imagesc(R,time,S2LC, [CFARthresh max(max(S2LC))]);colorbar
figure(34);
mesh(R(1:150),time(1:end-1),S2LC(:,1:150));
axis([0 R(150) 1 time(end) CFARthresh max(max(S2LC))])
title('RTI with smoothed mean range profile removed')
ylabel('time (s)');
xlabel('range (m)');
colorbar;

%% Plot range movie
if ~strcmp(rangeMovie,'no')
    Im = S2L;
    %     Im = S2;
    maxVal = max(max(Im));
    minVal = min(min(Im));
    [Rw,Cl]=size(Im);
    width = 500;
    delR = max_range/zpad;
    range = linspace(0,width*delR,width);
    figure(40)
    for ii=1:Rw
        plot(range, Im(ii,1:width))
        axis([0,width*delR,minVal maxVal])
        timeDisp = sprintf('time = %0.2f',time(ii));
        text(20,maxVal-10,timeDisp);
        time(ii);
        pause(0.1)
    end
end
%% Track targets
% if ~strcmp(trackTgt,'no')
% %     for row = 1:size(
% end
