function F0 = frequence(x,Fs,Fmin,Fmax,H)

% Assume already working on a given window
N = length(x);
% Fill in any missing arguments
if nargin < 5
    H = 4;
end
if nargin < 4
    Fmax = 900;
end
if nargin < 3
    Fmin = 50;
end

% Adapt Nfft to window size
Nfft = 2^nextpow2(N);
% Adjust spectral precision
dFfft = Fs/Nfft;
% Calculate DFT frequency vector
f = 0:dFfft:(Fs-(1/Nfft));
% Calculate DFT
X = fft(x,Nfft);
% Keep only useful part of DFT
X = X(1:(0.5*Nfft));
f = f(1:(0.5*Nfft));

% Initialise decimation matrix
XH = zeros(0.5*Nfft,H);
% Calculate progressive decimations
for k = 1:H
    XHk = abs(X(1:k:(0.5*Nfft)));
    R = length(XHk);
    XH(1:R,k) = XHk;
end

% Calculate spectral product
P = prod(XH,2);
% Keep only valid (non-zero) part of product
P = P(1:R);
% Adjust frequency vector to match
fP = f(1:R);
% Recalculate spectral precision
dFP = R/fP(R);
% Find indices corresponding to given Fmin and Fmax
Nmin = floor(dFP*Fmin);
Nmax = floor(dFP*Fmax);
% Search for max amplitude and corresponding index between Fmin and Fmax
if Nmax > length(P)
    Nmax = length(P);
    if Nmin >= Nmax
        Nmin = Nmax-1;
    end
end
[~,kF0] = max(P(Nmin:Nmax));
% Adjust frequency vector accordingly
fF0 = f(Nmin:Nmax);
% Find the fundamental (frequency value corresponding to max amplitude)
F0 = fF0(kF0);