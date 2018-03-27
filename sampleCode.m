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

for offset = 0:0.1:3
    
    %CFAR Implementation
    %2D window CA-CFAR
    refLength=32;
    guardLength=10;
    offset=1;
    cfarWin=ones((refLength+guardLength)*2+1,(refLength+guardLength)*2+1);
    cfarWin(refLength+1:refLength+1+2*guardLength,refLength+1:refLength+1+2*guardLength)=0;
    cfarWin=cfarWin/sum(cfarWin);

    figure
    for sampleNum = 1:16
        RangexDopplerxChannel_sum_over_channels = zeros(256,128);
        for ii = 1:4 %This is a loop over the 4 Rx channels
            %Pull out a single sample of data (for a single channel)
            FastFreqxSlowTime = squeeze(FastFreqxSlowTimexChannelxSample(:,:,ii,sampleNum)); %extract the current FastFreq and SlowTime for the current channel and sample

            RangexSlowTime =   fft(FastFreqxSlowTime.*Win2D,NFFT,1).*1/ScaWin; %Range by pulse
            RangexDopplerNonShifted  =   fft(RangexSlowTime.*WinVel2D, NFFTVel, 2)./ScaWinVel; %Range by Doppler
            RangexDopplerxChannel(:,:,ii) = fftshift(RangexDopplerNonShifted, 2); %Shift to properly align 0 Doppler to middle
            RangexDopplerxChannel_sum_over_channels = RangexDopplerxChannel_sum_over_channels + RangexDopplerxChannel(:,:,ii);

        end
            % Display range doppler map
            %Make image of data, note only half of the range is used (other half is invalid)


            noiseLevel=conv2(RangexDopplerxChannel_sum_over_channels(1:128,:),cfarWin,'same');
            cfarThreshold=noiseLevel+offset;

            RangexDopplerxChannel_sum_over_channels(1:128,:)= RangexDopplerxChannel_sum_over_channels(1:128,:) - cfarThreshold;

            imagesc(vVel, vRange(1:128), 10*log10(abs(RangexDopplerxChannel_sum_over_channels(1:128,:))));
            title('Range-doppler over all channels', num2str(offset))
            grid on;
            xlabel('v (m/s)');
            ylabel('R (m)');
            colormap('jet')
            caxis([-20 20])
            set(gca,'YDir','normal')
        pause(1.0)
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Search Angle Example                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Pull out data for a particular CPI (The 10th one in this case)
FastFreqxSlowTimexChannel = squeeze(FastFreqxSlowTimexChannelxSample(:,:,:,10));

RangexSlowTimexChannel = fft(FastFreqxSlowTimexChannel.*repmat(Win2D,1,1,4),NFFT,1).*1/ScaWin; %Range by pulse
RangexDopplerxChannel = fft(RangexSlowTimexChannel.*repmat(WinVel2D,1,1,4), NFFTVel, 2)./ScaWinVel; %Range by Doppler
%Shift the output of the fft so that 0 Doppler is in the middle (the 2
%indicates that the fftshift is done on the second dimension of the data (Doppler)
RangexDopplerxChannel = fftshift(RangexDopplerxChannel, 2);

figure
angleVec = linspace(-30, 30,61); %Vector of angles to search through
for ii = 1:numel(angleVec)
    c0 = 3e8;
    fc = 24.15e9;
    lambda = c0/fc;
    angle = angleVec(ii); %50;
    
    %Compute the phase for a wavefront coming in at the angle of interest
    %Note that there are 4 values and this assumes that the Rx elements are
    %spaced at 0.5 wavelengths (which should be close to true, might need
    %refined with actual values)
    p = (2*pi*sind(-angle)*[-1.5; -0.5; 0.5; 1.5]*(lambda/2))/lambda;
    %Compute the complex weights which will be applied to the 4 channels to
    %"steer" the array in the look direction and normalize the values
    w = exp(1i*p);
    w = w/norm(w);
    
    %Combine the four channels with the weights computed earlier
    RangexDoppler = abs(sum(RangexDopplerxChannel .*repmat(reshape(conj(w),1,1,4),nRng,nDopp,1),3));
    
    %Make image of data, note only half of the range is used (other half is invalid)
    imagesc(vVel, vRange(1:128), 10*log10(abs(RangexDoppler(1:128,:))));
    title(['Scan Angle = ',num2str(angle),' deg'])
    grid on;
    xlabel('v (m/s)');
    ylabel('R (m)');
    colormap('jet')
    caxis([-20 20])
    set(gca,'YDir','normal')
    pause(0.2)
end
