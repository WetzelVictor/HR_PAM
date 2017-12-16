%% Cosinus modulated bank filter
% Victor Wetzel
clear all; close all;
affich=1;
%% HARMONIC TEST SIGNAL
% const
Fe = 44100; Te = 1/Fe; N = Fe; N_FFT = 2^10;

% vectors
f = [200:400:20000]';
N_freq = length(f);
A = 2*rand(N_freq,1);

t = [0:N-1]'*Te;
DF = Fe/N_FFT;
FMAX_HZ = Fe/2;
F_VEC = [-FMAX_HZ+DF:DF:FMAX_HZ];

SIG_IN = zeros(N,1);

% signal
for i=1:N_freq
  SIG_IN = SIG_IN +A(i)*sin(2*pi*t*f(i));
end
SIG_IN = SIG_IN + 0.05*genAR(5,N)';
if affich
  figure
  plot(F_VEC, db(abs(fftshift(fft(SIG_IN,N_FFT)))))
  title('Spectre du signal d''entrée synthétique')
  xlabel('Fréquence (Hz)')
  ylabel('Amplitude (dB)')
  xlim([0, FMAX_HZ])
end

%% Synthesize a filter with remez algorithm
FMAX_HZ = Fe/2;
N_BANDES = 4;
CUTOFF_HZ = FMAX_HZ / N_BANDES;
CUTOFF_BINS = CUTOFF_HZ / (Fe*2); % Fréquence réduite
%Pourquoi Fe*2? Sinon les filtres se chevauchent

% Conception du filtre
d1 = 10e-6;
d2 = 10e-5;

DELTA = 0.9;
nuC = CUTOFF_BINS * (1 - DELTA);
nuA = CUTOFF_BINS * (1 + DELTA);

L = floor(2/(3*(nuA - nuC)) * log10(1/(10*d1*d2))); % Ordre
 % calcul de la réponse impulsionnelle
h = firpm(L - 1, 2*[0 nuC nuA 0.5], [1 1 0 0], [d2 d1]);

%% SYNTHESIZE FILTER BANK
[SIG_BANKED,Vm,H] = CM_filterbank(h,h,N_BANDES,SIG_IN);

if affich
  figure
  plot(F_VEC, db(abs(fftshift(fft(H',N_FFT)))))
  hold on
  plot(F_VEC, db(sum(abs(fftshift(fft(H',N_FFT))),2)))
  title('Banc de filtre')
  xlabel('Fréquence (Hz)')
  ylabel('Amplitude (dB)')
  legend('Filtre 1','Filtre 2','Filtre 3','Filtre 4', 'Somme')
  xlim([0, FMAX_HZ])
end

%% FILTERING SIGNAL
X = [];
X = repmat(zeros(length(conv(H(1,:)',SIG_IN)),1),1,N_BANDES);

for i=1:N_BANDES
  X(:,i) = conv(H(i,:).',SIG_IN)';
end

if affich
  figure
  plot(F_VEC, db(abs(fftshift(fft(X,N_FFT,1)))))
  hold on
  plot(F_VEC, db(abs(fftshift(fft(SIG_IN,N_FFT)))))
  title('Banc de filtre')
  xlabel('Fréquence (Hz)')
  ylabel('Amplitude (dB)')
  legend('Filtre 1','Filtre 2','Filtre 3','Filtre 4', 'Original')
  xlim([0, FMAX_HZ])
end
