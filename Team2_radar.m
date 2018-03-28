% add paths to the provided functions/classes (note, this is really only
% used in this example for some of the weighting functions and for pulling
% the sampling rate and number of samples per ramp)

addpath('./Support/DemoRadUsb');
addpath('./Support/Class');

%Load in the pre-recorded data
load('OutboundCarSample.mat')

% Calculate range/Doppler
Brd         =   Adf24Tx2Rx4();
N           =   Brd.Get('N'); %256
fs          =   Brd.Get('fs'); %1000000
c0          =   3e8; %speed of light

NFFT = 256;
NFFTVel = 128;

% Processing of range profile and setting up FFT weights
Win2D           =   Brd.hanning(N,Cfg.Np); %range taper
ScaWin          =   sum(Win2D(:,1)); %Value used to normalize the fft

kf              =   (Cfg.fStop - Cfg.fStrt)/Cfg.TRampUp; %slope of chirp
vRange          =   (0:NFFT-1).'./NFFT.*fs.*c0/(2.*kf); %range values after fft (meters)
fc              =   (Cfg.fStop + Cfg.fStrt)/2; %center frequency

WinVel          =   Brd.hanning(Cfg.Np); %Velocity taper
ScaWinVel       =   sum(WinVel); %Value used to normalize the fft
WinVel2D        =   repmat(WinVel.',numel(vRange),1);

%Compute Doppler frequencies for values after fft is computed
vFreqVel        =   (-NFFTVel/2:NFFTVel/2-1).'./NFFTVel.*(1/Cfg.Tp);
%Convert Doppler frequencies into velocity values
vVel            =   vFreqVel*c0/(2.*fc); 

%loop over all 4 channels
nRng = 256;
nDopp = 128;
%Initialize matrix for holding range/Doppler maps for each channel
RangexDopplerxChannel = zeros(nRng,nDopp,4);

% Dont touch above code, initializes data sets and computes the velocities
% for us


    
%CFAR Implementation
%2D window CA-CFAR
%https://www.mathworks.com/matlabcentral/answers/165561-how-to-write-a-m-file-code-to-cfar-for-fmcw-radar
refLength=10;
guardLength=2;
offset=0.05;

cfarWin=ones((refLength+guardLength)*2+1,(refLength+guardLength)*2+1); %2D window

%need to modify kernel for 2D window
cfarWin(refLength+1:refLength+1+2*guardLength,refLength+1:refLength+1+2*guardLength)=0;

cfarWin=cfarWin./sum(cfarWin(:));


% belief is used to track the object
belief = ones(128, 128) / (128*128);

% the range kernel defines how much noise we expect in our motion model
sigma = 3;
sz = 30;    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
range_kernel = exp(-x .^ 2 / (2 * sigma ^ 2))';
range_kernel = range_kernel / sum (range_kernel); % normalize

% the doppler kernel defines how much acceleration we can expect
sigma = 30;
sz = 100;    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
doppler_kernel = exp(-x .^ 2 / (2 * sigma ^ 2))';
doppler_kernel = doppler_kernel / sum (doppler_kernel); % normalize

dt = 1;
figure(1);
for sampleNum = 1:16
    accum = zeros(256,128);
    for ii = 1:4 %This is a loop over the 4 Rx channels
        %Pull out a single sample of data (for a single channel)
        FastFreqxSlowTime = squeeze(FastFreqxSlowTimexChannelxSample(:,:,ii,sampleNum)); %extract the current FastFreq and SlowTime for the current channel and sample

        RangexSlowTime =   fft(FastFreqxSlowTime.*Win2D,NFFT,1).*1/ScaWin; %Range by pulse
        RangexDopplerNonShifted  =   fft(RangexSlowTime.*WinVel2D, NFFTVel, 2)./ScaWinVel; %Range by Doppler
        RangexDopplerxChannel(:,:,ii) = fftshift(RangexDopplerNonShifted, 2); %Shift to properly align 0 Doppler to middle
        accum = accum + RangexDopplerxChannel(:,:,ii);
    end
    
    % Calculate range doppler map
    % Make image of data, note only half of the range is used (other half is invalid)
    pmf = abs(accum(1:128,:)/4);
    noiseLevel=conv2(pmf,cfarWin,'same');
    cfarThreshold=noiseLevel+offset;

    % create binary mask from cfarTheshold
    Idx = pmf - cfarThreshold <= 0;
    % remove detections around 0 velocity
    Idx(:,vVel < 0.5 & vVel > - 0.5) = 1;
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
    
    % Tracking algorithm
    
    % Prediction update - applies motion model
    for i = 1:size(belief, 2)
        % shift ranges according to velocity
        if vVel(i) > 0
            n = floor(vVel(i) * dt + 0.5);
            belief(:, i) = [zeros(n, 1); belief(1:end-n, i)];     
        else
            n = floor(-vVel(i) * dt + 0.5);
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
    
    % Display range doppler map
    subplot(2, 2, 1);
    imagesc(vVel, vRange(1:128), 10*log10(abs(accum(1:128,:)/4)));
    title('Raw Range-doppler over all channels')
    grid on;
    xlabel('v (m/s)');
    ylabel('R (m)');
    colormap('jet')
    colorbar;
    %caxis([-20 20])
    set(gca,'YDir','normal')
    
    subplot(2, 2, 2);
    imagesc(vVel, vRange(1:128), pmf);
    title('Range-doppler CFAR')
    grid on;
    xlabel('v (m/s)');
    ylabel('R (m)');
    colormap('jet')
    colorbar;
    caxis([0 1])
    set(gca,'YDir','normal')
    
    subplot(2, 2, 3);
    imagesc(vVel, vRange(1:128), belief);
    title('Histogram filter');
    colorbar;
    set(gca,'YDir','normal')
    
    [i, j] = find(belief == max(belief(:)));
    subplot(2, 2, 4);
    plot(vVel(j), vRange(i), 'r*');
    xlim([vVel(1) vVel(end)]);
    ylim([vRange(1) vRange(128)]);
    title('Tracker');
    colorbar;
    pause(1.0)
end

