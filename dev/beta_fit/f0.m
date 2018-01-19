function [Xk, freq0] = f0(x, Fs, dF, Fmin, Fmax, H)
% Estimation de F0 par la méthode du produit spectral
%
% INPUT:
%    - x: signal d'entrée
%    - Fs: FrÃ©quence d'échnatillonage
%    - dF: pas en fréquence? (def =Fs/N)
%    - Fmin/Fmax: Bornes de recherche (def: [50, 900])
%    - H: ordre de contraction du spectre (def: H = 4)
%
% OUTPUT: 
%    - F0: frÃ©quence fondamentale estimÃ©e (Hz)

%% DEFAULT ARGUMENTS

if nargin < 6
  H = 4;
end
if nargin < 5
  Fmax = 900;
end
if nargin < 4
  Fmin = 50;
end
if nargin < 3
  dF = Fs/length(x);
end

%% INIT
N = length(x);
Nfft = Fs/dF; % ordre de la TFD

% Calcul ordre en fonction de la precision dF
% p = nextpow2(Nfft);
% Nfft=2^p;
fmax = Fs/2;

% Valeur maximale de R
R = round(Nfft/(2*H) + 1);
fRmax = Fs/Nfft*R;
f = [1/R:1/R:1]*fRmax;

% FFT
Xk = fft(x,Nfft);

%% Matrice Transposition
transMat = zeros(H, R);
transMat(1,:) = abs(Xk(1:R));

for i=2:H
  for j=1:R
        transMat(i,j) = abs(Xk(1+H*(j-1)));
  end
end

%% Produit spectral
specProd = prod(transMat,1);

% Recherche du maximum
[maxP, indP] = max(specProd);

% Fréquence
freq0 = f(indP);
freq0 = round(freq0);

%% Return conditions
if freq0 < Fmin
    freq0 = Fmin;
elseif freq0 > Fmax
    freq0 = Fmax;
end
end
