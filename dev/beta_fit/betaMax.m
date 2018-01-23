clear all, close all
[x,Fs] = audioread('sinus_45deg_P2.wav'); % Load audio data

N = floor(0.2*Fs); % Number of samples in 0.25 seconds (analysis length)
Nstart = floor(0.2*Fs); % Start analysing at 0.1 seconds (cut attack)
dF = Fs/N; % Discrete frequency step

% Truncate and window audio data for analysis
x_fen = x(Nstart:(N+Nstart-1)).*hann(N);

Nfft = 2^nextpow2(N); % Determine DFT points from analysis length
dFfft = Fs/Nfft; % DFT frequency step

X = fft(x_fen,Nfft); % Spectrum of truncated + windowed audio data
X = X(1:(0.5*Nfft)); % Take real part only

F0 = frequence(x_fen,Fs); % Determine fundamental frequency
W = db(fft(hann(N))); % DFT of Hann window in dB
nw = find((W<=(max(W)-60)),1); % Index of first point 60 dB below max
Ww = 2*nw*dF; % Hann window principal lobe width in frequency
kF0 = (F0/dFfft)+1; % Index of fundamental

betaVect = 0:0.0001:0.1-0.0001; % Vector of possible beta values
alpha = 0.01; % Determined empirically
E = zeros(length(betaVect),1); % Vector to store total peak energies
fp = [];
for b = 1:length(betaVect)
    beta = betaVect(b);
    Xsup = abs(X); % Initialise amplitude spectrum to suppress
    XF0 = Xsup(kF0); % Amplitude of fundamental frequency
    k1 = round(((F0-Ww)/dFfft)+1); % Index of lower bound
    k2 = round(((F0+Ww)/dFfft)+1); % Index of upper bound
    tF0 = min(Xsup(k1:k2)); % Threshold to use as while loop condition
    Xfl = XF0; % Initialise current peak amplitude
    Xsup(k1:k2) = tF0; % Suppress fundamental
    E(b) = Xfl; % Initialise total peak energy
    l = 2; % Next partial to analyse (first in loop)
    while Xfl > tF0
        finharmo = l*F0*sqrt(1+((l^2)-1)*beta); % Theoretical value
        flmin = (1-alpha)*finharmo; % Lower bound frequency
        flmax = (1+alpha)*finharmo; % Upper bound frequency
        fl = frequence(x_fen,Fs,flmin,flmax); % "Fundamental" of partial            
        k1 = round(((fl-Ww)/dFfft)+1); % Lower bound index
        k2 = round(((fl+Ww)/dFfft)+1); % Upper bound index
        Xfl = max(Xsup(k1:k2)); % Amplitude of peak
        Xsup(k1:k2) = min(Xsup(k1:k2)); % Suppress partial
        E(b) = E(b)+Xfl; % Add to total peak energy
        l = l+1; % Move to next partial
    end
    fp = [fp fl];
end

[~,e] = max(E); % Index of max peak energy
maxBeta = betaVect(e); % Optimized beta value