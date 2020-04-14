%=============================================
%SMS-Matlab like emulation
clear all
close all
%==== USER DATA =====
DAFx_in = audioread('273175__xserra__la-vaca-cega-xavier.wav'); % wave file
SR = 44100; %Ssampling rate
w1Length = 2048; % analysis window size
n1 = 256; % analysis window hop size
nPeaks = 100; % number of peaks detected
nSines = 50; % number of sinuosoids to track (and synthetise)
minSpacePeaks = 2; % minimum space (bins) between two picked peaks
rgain = 1.; % gain for the residual component
MaxFreq = 11000; % maximum frequency, in Hertz, for plottings
MinMag = -100; % minimum magnitude, idnB , for plottings
zp = 2; % zero-padding coefficient
%--- figure data
%fig1 = ’yes’; % if uncommented, will plot the Blackman-Harris
%fig2 = ’yes’; % if uncommented, will plot the peaks detection
%fig3 = ’yes’; % if uncommented, will plot the peak trackings
%fig4 = ’yes’; % if uncommented, will plot the original and
%fig5 = ’yes’; % if uncommented, will plot the peak trackings
% only at the eonfd the process
%fig6 = ’yes’; % if uncommented, will plot the original signal,
% its sine and residual part,
% and the transformed signal
% window
% and tracking in one frame
% real-time
% the transformedF FT in one frame
%=== Definition of the Windows===

%--- definition of the analysis window
fConst=2*pi/(w1Length+1-1);
w1 = [1:w1Length]';
w1=.35875 -.48829*cos(fConst*w1)+.14128*cos(fConst*2*w1)- .01168*cos(fConst*3*w1);
w1= w1/sum(w1)*2;
N= w1Length*zp; %new size of the window
%--- synthesis window
w2 = w1;
n2 = n1;
%--- triangular window
wt2=triang(n2*2+1); % triangular window
%--- main lobe table of bh92
[bh92SINE2SINE,bh92SINE2SINEsize]=bh92SINE2SINEgeneration;
%--- data for the loops
frametime = n1/SR;
pin = 0;
pout= 0;
TuneLength=length(DAFx_in);
pend=TuneLength-w1Length;

%=== Definition of the data arrays ===
DAFx_in = [zeros(w1Length/2-n1-1 ,2); DAFx_in];
DAFx_outsine = zeros(TuneLength,1);
DAFx_outres = zeros(TuneLength,1);
iloc = zeros(nSines, 1);
ival = zeros(nSines, 1);
iphase = zeros(nSines,1);
previousiloc = zeros(nSines, 1);
previousival = zeros (nSines, 1);
maxSines = 400; % maximum voices for harmonizer
syniloc = zeros (maxSines, 1);
synival = zeros(maxSines, 1);
previoussyniloc = zeros(maxSines,1);
previousiphase = zeros(maxSines,1);
currentiphase = zeros(maxSines, 1);
%--- arrays for the sinus' frequencies and amplitudes
SineFreq = zeros(nSines,ceil(TuneLength/n2));
SineAmp = zeros(nSines,ceil(TuneLength/n2));
pitch = zeros(1,1+ceil(pend/n1));
pitcherr = zeros(1,1+ceil(pend/n1));
%--- arrays for the partial tracking
%--- arrays for the sinus) frequencies and amplitudes
%--- creating figures ---
if (exist('fig1'))
    h = figure(1); set (h, 'position', [10, 45, 200, 200]);
end
if (exist('fig2'))
    h = figure(2) ; set(h, 'position', [10, 320, 450, 350]);
    axisFig2 = [0 MaxFreq MinMag 0]; zoom on;
end
if(exist('fig3'))
    h = figure(3); set(h, 'position', [220, 45, 550, 200]);
    axisFig3 = [l l+ceil(pend/n1) 0 MaxFreq]; zoom on;
end
if(exist('fig4'))
    h = figure(4); set(h, 'position', [470, 320, 450, 350]);
    axisFig4 = [0 MaxFreq MinMag 0]; zoom on;
end

if(exist('fig5'))
    h = figure(5) ; set(h, 'position', [220, 45, 550, 200]);
    axisFig5 = [l l+ceil(pend/n1) 0 MaxFreq] ; zoom on;
end

%--- plot the Blackman-Harris window
if (exist('fig1'))
    figure (1)
    plot(20*log10(abs(fftshift(fft(bh92SINE2SINE)/bh92SINE2SINEsize))))
    title('B1ackman-Harris window'); xlabel('Samples');
    ylabel('Amplitude')
end

tic
%UUUUUUUUUUUUUUUUUUUUUUUUUU
disp('analyzing frame ...') ;

while pin<pend
    %--- windowing
    grain = DAFx_in(pin+1:pin+w1Length).*w1(1:w1Length);
    %--- zero padding
    padgrain = zeros(N, 1);
    padgrain(1:w1Length/2) = grain(w1Length/2+1:w1Length);
    padgrain(N-w1Length/2+1 :N) = grain(1: w1Length/2);
    %--- fft computation
    f = fft(padgrain);
    r = abs(f);
    phi = angle(f);
    ft = r.*exp(j*phi);
    
    %===== Analysis =====
    
    %--- peak detection (and their plottings)
    [ftloc, ftval]=pickpeaks(r(1:N/2),nPeaks,minSpacePeaks);
    %--- calculate interpolated values (peak position, phase, amplitude)
    [iftloc, iftphase, iftval]= interpolatedvalues(r, phi, N, zp, ftloc, ftval);
    
    %--- pitch detection
    [pitchvalue, pitcherror, isHarm]= pitchDetection(r, N, SR, nPeaks, iftloc, iftval);
    pitch(1+pin/n1)=pitchvalue*isHarm;
    pitcherr(1+pin/n1)= pitcherror;
    
    %--- peaks tracking
    if(pin==0) %--- for the first frame
        nNewPeaks = nSines;
    else %--- creating new born tracks
        for i=1:nSines
            if (previousiloc(i)==0)
                [previousiloc(i), previousival(i)]= CreateNewTrack(iftloc, iftval, previousiloc, previousival, nSines, MinMag);
                nNewPeaks= nNewPeaks-1;
            end
        end
        
        %--- simple Peak tracker
        [iloc, ival, iphase, previousiloc, previousival, distminidex]= peakTrackSimple(nSines, nPeaks, N, SR, pitchvalue, iftloc, iftval, iftphase, isHarm, previousiloc, previousival);
    end
    
    %--- savings
    previousival= ival;
    previousiloc= iloc;
    SineFreq(:, 1+pin/n1)=max((iloc-1)/N*SR, 0.);
    %frequency of the partials
    SineAmp(:,1+pin/n1)=max(ival, MinMag);
    %amplitudes of the partials
    syniloc(1:nSines)=max(1, iloc);
    synival(1:nSines)= ival;
    
    if(exist('fig3')) %plot: the trackings of partials
        figure(3); clf; hold on
        PlotTracking(SineFreq(:,1:1+pin/n1), pitch(1:1+pin/n1));
        xlabel('Frame Number'); ylabel('Frequency (Hz)');
        axis(axisFig3);title('Peak tracking'); drawnow
    end
    
    %---residual computation
    resfft = ft;
    if(isHarm==1)
        resfft= resfft-sinefillspectrum(iloc, ival, iphase, nSines, w1length, zp, bh92SINE2SINE, bh92SINE2SINEsize);
    end
    
    %--- figures
    if(exist('fig2'))
        figure(2); clf; hold on
        %plot: FFT of the windowed signal (Hz, dB)
        plot((1:N/2)/N*SR, 20*log10(r(1:N/2)));
        for l=1:nPeaks    %plot: the peaks detected
            plot([ftloc(1)-1 ftloc(1)-1]/N*SR, [20*log10(ftval(1)), MinMag-1], 'r:x');
        end
        for l=1:nSines   %plot: sines tracked and the residual part
            plot([iloc(1)-1, iloc(1)-1]/N*SR, [ival(1), MinMag-1], 'k')
        end
        plot([iloc(1)-1, iloc(1)-1]/N*SR, [ival(1), MinMag-1], 'k')
        if(isHarm)   %plot: true pitch of each harmonic
            for l=1:nSines
                plot([pitchvalue*1, pitchvalue*1], [1, MinMag-1], 'y')
            end
        end
        xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); axis(axisFig2); title('Peak detection and tracking for one frame'); drawnow
    end
    
    nSynSines= nSines;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %=====Transformations =====
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %===== gender change: woman to man =====
    tr='m2f'; %male to female
    %tr='f2m' ; %female to male
    if (isHarm == I)
        pitchmin=lOO;
        pitchmax=500;
        sssmax = 50;
        if (pitchvalue<pitchmin)
            sss = 0;
        elseif (pitchvalue>pitchmax)
            sss = sssmax;
        else
            sss = (pitchvalue-pitchmin)/((pitchmax-pitchmin)/sssmax);
        end
        if (tr=='f2m')
            sss = -sss;
            pt = .5;
        else
            pt= 2 ;
        end
        %--- spectral shape computation
        [spectralShape,shapePos]=CalculateSpectralShape(iloc,ival,MinMag,N);
        %--- spectral shape shift
        syniloc = zeros(nSines,l);
        if shapePos > 1
            [shiftedSpectralShape,shapePos] = SpectralShapeShift(sss,iloc, ival, spectralShape, shapePos, N, SR);
        end
        syniloc = iloc;
        %linear interpol. of the spectral shape for synival computation
        if shapePos > 1
            synival = interpl(shiftedSpectralShape(l,l:shapePos+l)', shiftedSpectralShape(2,1:shapePos+l)', syniloc, 'linear');
        else
            synival = ival;
        end
        %--- pitch transposition
        pt = 0.5;
        [syniloc, synival] = PitchTransposition(iloc,ival,spectralShape, shapePos,pt ,N) ;
        %--- comb filtering the residual
        CombCoef = 1;
        if (isHarm==l)
            resft = CombFilter(resft, N , SR/( pitchvalue*pt) , CombCoef) ;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %===== Synthesis =====
    if(pin>0)
        for i=1: nSynSines
            if (syniloc(i)~=0)
                ifreq= (previoussyniloc(distminindex(i)+syniloc(i))/2);
                %average bin
                freq= (ifreq-1)/N*SR; %freq in Hz (if loc=1 --> freq=0)
                currentiphase(i)=unwrap2pi(previousiphase(distminindex(i))+2*pi*freq*frametime);
            end
        end
    end
    
    previoussynival= synival;
    previoussyniloc= syniloc;
    previousiphase= currentiphase;
    
    %--- compute sine spectrum
    padsynthft= sinefillspectrum(syniloc, synival, currentiphase, nSynSines, w1Length, zp, bh92SINE2SINE, bh92SINE2SINEsize);
    if(isHarm==0)
        padsynthft= zeros(size(padsynthft));
    end
    
    %---residual computation
    respadgrain= real(ifft(resfft));
    resgrain= [respadgrain(N-w1Length/2+1:N); respadgrain(1:w1Length/2)]./w2(1:w1Length);
    ressynthgrain=wt2(1:n2*2).*resgrain(w1Length/2-n2:w1Length/2+n2-1);
    DAFx_outres(pout+1:pout+n2*2)= DAFx_outres(pout+1:pout+n2*2)+resynthgrain;
    
    %--- sinusoidal computation
    sinpadgrain= real(ifft(padsynthft));
    singrain=[sinpadgrain(N-w1Length/2+1:N); sinpadgrain(1:w1Length/2)]./w2(1:w1Length);
    sinsynthgrain= wt2(1:n2*2).*singrain(w1Length/2-n2:w1Length/2+n2-1);
    DAFx_outsine(pout+1:pout+n2*2)= DAFx_outsine(pout+1:pout+n2*2)+sinsynthgrain;
    
    %--- figure with original signal and transformed signal FFT
    synthr=abs(fft(respadgrain+sinpadgrain));
    if(exist('fig4'))
        figure(4); clf; hold on
        plot((1:N/2)/N*SR, 20*log10(r(1:N/2)), 'b:'); axis(axisFig4);
        plot((1:N/2)/N*SR, 20*log10(synthr(1:N/2)), 'r');
        figure(4);
        xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); axis(axisFig4);
        title('FFT of original (blue) and transformed (red) signals');
        drawnow
    end
    
    %--- increment loop indexes
    pin= pin+n1;
    pout= pout+n2;
    disp(pin/n1);
end
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUU
toc

%===== write output sounds =====
DAFx_in= DAFx_in(w1Length/2-n1:length(DAFx_in));
%remove the zeros added for the process
DAFx_outresynth= DAFx_outsine(1:TuneLength)+rgain*DAFx_outres(1:TuneLength);
mm=max(abs(DAFx_outresynth));
audiowrite(DAFx_outresynth/mm, SR, 'DAFx_out.wav');
audiowrite(DAFx_outsine/mm, SR, 'DAFx_outsine.wav');
audiowrite(DAFx_outres/mm, SR, 'DAFx_outres.wav');

if(exist('fig3')==0 & exist('fig5')) %plot:trackings of partials
    %only at the end of the process
    figure(5); clf; hold on
    PlotTracking(SineFreq(:,1:1+pend/n1), pitch(1:1+pend/n1));
    xlabel('Frame number'); ylabel('Frequency (Hz)'); axis(axisFig5);
    title('Peak tracking'); drawnow
end

if(exist('fig6')) %plot the input signal, its sinus
    % and its residual part, and the transformed
    figure(6)
    subplot(4,1,1); plot(DAFx_in); xlabel('input signal');
    subplot(4,1,2); plot(DAFx_outsine); xlabel('sinus part');
    subplot(4,1,3); plot(DAFx_outres); xlabel('residual part');
    subplot(4,1,4); plot(DAFx_outresynth); xlabel('resynthesized signal');
end

    
