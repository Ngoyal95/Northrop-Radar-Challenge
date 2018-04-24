% AN24_06 -- Range-Doppler basics
% Evaluate the range-Doppler map for a single channel
close all;

% (1) Connect to DemoRad: Check if Brd exists: Problem with USB driver
% (3) Configure RX
% (4) Configure TX
% (5) Start Measurements
% (6) Configure calculation of range profile and range doppler map for
% channel 1

c0          =   3e8;
%--------------------------------------------------------------------------
% Include all necessary directories
%--------------------------------------------------------------------------
CurPath = pwd();
addpath([CurPath,'/../../DemoRadUsb']);
addpath([CurPath,'/../../Class']);

%--------------------------------------------------------------------------
% Setup Connection
%--------------------------------------------------------------------------
Brd         =   Adf24Tx2Rx4();

Brd.BrdRst();

%--------------------------------------------------------------------------
% Software Version
%--------------------------------------------------------------------------
Brd.BrdDispSwVers();

%--------------------------------------------------------------------------
% Configure Receiver
%--------------------------------------------------------------------------
Brd.RfRxEna();
TxPwr           =   80;
NrFrms          =   4;

%--------------------------------------------------------------------------
% Configure Transmitter (Antenna 0 - 4, Pwr 0 - 31)
%--------------------------------------------------------------------------
Brd.RfTxEna(1, TxPwr);

%--------------------------------------------------------------------------
% Configure Up-Chirp
% TRamp is only used for chirp rate calculation!!
% TRamp, N, Tp can not be altered in the current framework
% Effective bandwidth is reduced as only 256 us are sampled from the 280 us
% upchirp.
% Only the bandwidth can be altered
% 
% The maximum number of chirps is StrtIdx = 0 and StopIdx = 128
%--------------------------------------------------------------------------
Cfg.fStrt       =   24.0e9;
Cfg.fStop       =   24.3e9;
Cfg.TRampUp     =   280e-6;
Cfg.Tp          =   284e-6;
Cfg.N           =   256;
Cfg.StrtIdx     =   0;
Cfg.StopIdx     =   64;
Cfg.Np          =   Cfg.StopIdx - Cfg.StrtIdx;

Brd.RfMeas('Adi', Cfg);

%--------------------------------------------------------------------------
% Read actual configuration
%--------------------------------------------------------------------------
NrChn           =   Brd.Get('NrChn');
N               =   Brd.Get('N');
fs              =   Brd.Get('fs');

%--------------------------------------------------------------------------
% Configure Signal Processing
%--------------------------------------------------------------------------
% Processing of range profile
Win2D           =   Brd.hanning(N,Cfg.Np);
ScaWin          =   sum(Win2D(:,1));
NFFT            =   2^12;
NFFTVel         =   2^8;
kf              =   (Cfg.fStop - Cfg.fStrt)/Cfg.TRampUp;
vRange          =   [0:NFFT-1].'./NFFT.*fs.*c0/(2.*kf);
fc              =   (Cfg.fStop + Cfg.fStrt)/2;

RMin            =   0.5;
RMax            =   10;

[Val RMinIdx]   =   min(abs(vRange - RMin));
[Val RMaxIdx]   =   min(abs(vRange - RMax));
vRangeExt       =   vRange(RMinIdx:RMaxIdx);

WinVel          =   Brd.hanning(Cfg.Np);
ScaWinVel       =   sum(WinVel);
WinVel2D        =   repmat(WinVel.',numel(vRangeExt),1);

vFreqVel        =   [-NFFTVel./2:NFFTVel./2-1].'./NFFTVel.*(1/Cfg.Tp);
vVel            =   vFreqVel*c0/(2.*fc); 

    
%CFAR Implementation
%2D window CA-CFAR
%https://www.mathworks.com/matlabcentral/answers/165561-how-to-write-a-m-file-code-to-cfar-for-fmcw-radar
refLength=10;
guardLength=2;
offset=0.000005;

cfarWin=ones((refLength+guardLength)*2+1,(refLength+guardLength)*2+1); %2D window
cfarWin(refLength+1:refLength+1+2*guardLength,refLength+1:refLength+1+2*guardLength)=0;

cfarWin=cfarWin./sum(cfarWin(:));


% belief is used to track the object
belief = ones(279, 256) / (279*256);

% the range kernel defines how much noise we expect in our motion model
sigma = 3;
sz = 30;    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
range_kernel = exp(-x .^ 2 / (2 * sigma ^ 2))';
range_kernel = range_kernel / sum (range_kernel); % normalize

% the doppler kernel defines how much acceleration we can expect
sigma = 60;
sz = 200;    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
doppler_kernel = exp(-x .^ 2 / (2 * sigma ^ 2))';
doppler_kernel = doppler_kernel / sum(doppler_kernel); % normalize

dt = 0;

%NEED TO CHANGE THIS TO INF. LOOP?
tic
while true
    accum = zeros(279,256);
%     datpath='RD/';
%     string=[datpath 'RD_' num2str(iter) '.mat'];
%     load(string)
%     
%     Datan        =   sum(Data,2);   
%     
%     MeasChn     =   reshape(Datan,N,Cfg.Np);
%     % Calculate range profile including calibration
%     RP          =   fft(MeasChn,NFFT,1);%/ScaWin;
%     RPExt       =   RP(RMinIdx:RMaxIdx,:);    

    Data        =   Brd.BrdGetData();     
    dt = toc;
    MeasChn     =   reshape(Data(:,1),N,Cfg.Np);
    % Calculate range profile including calibration
    RP          =   fft(MeasChn.*Win2D,NFFT,1).*Brd.FuSca/ScaWin;
    RPExt       =   RP(RMinIdx:RMaxIdx,:);    

    RD          =   fft(RPExt.*WinVel2D, NFFTVel, 2)./ScaWinVel;
    RD          =   fftshift(RD, 2);
    
    % Make image of data, note only half of the range is used (other half is invalid)
    pmf = abs(RD);
    noiseLevel=conv2(pmf,cfarWin,'same');
    offset = 1e-6;
    cfarThreshold=noiseLevel+offset;
    Idx = pmf - cfarThreshold <= 0;
    % remove detections around 0 velocity
    Idx(:,vVel < 1 & vVel > -1) = 1;
    % dilate the mask to get more information about detections
    Idx = logical(1 - Idx);
    se = strel('disk', 5);
    Idx = imdilate(Idx, se);
    Idx = logical(1 - Idx);

    % normalize data and apply mask
    pmf = 10*log10(pmf);
    pmf = pmf - min(pmf(:));
    pmf = pmf ./ max(pmf(:));
    pmf(Idx) = 0.001;
%     
%     % Tracking algorithm
    
    % Prediction update - applies motion model
    for i = 1:size(belief, 2)
        % shift ranges according to velocity
        if vVel(i) > 0
            n = min(floor(vVel(i) * dt + 0.5), 278);
            belief(:, i) = [zeros(n, 1); belief(1:end-n, i)];     
        else
            n = min(floor(-vVel(i) * dt + 0.5), 278);
            belief(:, i) = [belief(n+1:end, i); zeros(n, 1)];
        end
        
        % convolve with range_kernel to model uncertainty in position
        belief(:, i) = conv(belief(:, i), range_kernel, 'same');
      
    end
    
    % convolve with doppler_kernel to model uncertainty in velocity
    for i = 1:size(belief, 1)
        belief(i, :) = conv(belief(i, :), doppler_kernel, 'same');
    end
    
    
    % Correction update - uses range doppler map 
    for i = 1:size(belief, 1)
        for j = 1:size(belief, 2)
            belief(i, j) = belief(i, j) * pmf(i, j);
        end
    end
    belief = belief / sum(belief(:));
   
%     subplot(2, 2, 1);
%     %imagesc(vVel, vRange(1:128), 10*log10(abs(accum(1:128,:)/4)));
%     imagesc(vVel, vRangeExt, abs(RD));
%     title('Raw Range-doppler over all channels')
%     grid on;
%     xlabel('v (m/s)');
%     ylabel('R (m)');
%     colormap('jet')
%     colorbar;
%     %caxis([-20 20])
%     set(gca,'YDir','normal')
    
    %subplot(2, 2, 2);
    %imagesc(vVel, vRangeExt, pmf);
    %title('Range-doppler CFAR')
    %grid on;
    %%xlabel('v (m/s)');
    %ylabel('R (m)');
%     colormap('jet')
%     colorbar;
%     caxis([0 1])
%     set(gca,'YDir','normal')
%     
    %subplot(2, 2, 3);
    %imagesc(vVel, vRangeExt, belief);
    %title('Histogram filter');
    %colorbar;
    %set(gca,'YDir','normal')
    
    [i, j] = find(belief == max(belief(:)));
    
    %subplot(2, 2, 4);
    %plot(vVel(j), vRange(i), 'r*'); %plot the tracker (v,R)
    %xlim([vVel(1) vVel(end)]);
    %ylim([vRange(1) vRangeExt(end)]);
    %title('Tracker');
    %colorbar;
    if length(i) == 1
        fprintf('write status %d', write_detections(vRange(i), vVel(j)));
        fprintf('TRACKER: R: %f\t v: %f\n', vRange(i), vVel(j));
    end
    %pause(0.001)
end

