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

for offset = 1
    
    %CFAR Implementation
    %2D window CA-CFAR
    %https://www.mathworks.com/matlabcentral/answers/165561-how-to-write-a-m-file-code-to-cfar-for-fmcw-radar
    refLength=10;
    guardLength=2;
    offset=0.2;
    
    cfarWin=ones((refLength+guardLength)*2+1,(refLength+guardLength)*2+1); %2D window
    
    %need to modify kernel for 2D window
    cfarWin(refLength+1:refLength+1+2*guardLength,refLength+1:refLength+1+2*guardLength)=0;
   
    cfarWin=cfarWin./sum(cfarWin(:));

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

    
            pmf = abs(RangexDopplerxChannel_sum_over_channels(1:128,:));
            noiseLevel=conv2(pmf,cfarWin,'same');
            cfarThreshold=noiseLevel+offset;
            
            Idx = pmf - cfarThreshold <= 0;
            
            pmf = 10*log10(pmf);
            
            pmf = pmf - min(pmf(:));
            pmf = pmf ./ max(pmf(:));
            pmf(Idx) = 0;
            pmf(:,vVel < 0.5 & vVel > - 0.5) = 0;
            

            %pmf(pmf - cfarThreshold <= 0) = 0;
            
           
            %pmf = pmf * 100;
            imagesc(vVel, vRange(1:128), pmf);
            %title('Range-doppler over all channels', num2str(offset))
            grid on;
            xlabel('v (m/s)');
            ylabel('R (m)');
            colormap('jet')
            colorbar;
            caxis([0 1])
            set(gca,'YDir','normal')
        pause(1.0)
    end
end
