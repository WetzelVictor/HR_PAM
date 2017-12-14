%% Dev blanchiement
% Victor Wetzel - 2017
clear all; close all;

%% Initialisation
Fe = 44100;
N = 4096;

t = [0:N-1]/Fe;


sig = rand(N,1);
N_fft = 2^15;
win = hann(N);
N_median = 2^12;


fmax = Fe/2;
df = Fe/N_fft;
f_vec = [-fmax:df:fmax-df];

f = [250 500 1000 2000 4000 8000];
%% Signaux
% Signal de synthèse
a = [1.125 0.986 0.12345 0.2345];
b = [0.2 0.4];

for i=1:length(f)    
    sig = sig + sin(2*pi*t*f(i))';
end

sig = filter(1,a,sig);

% fft du signal test
fft_sig = fft(win.*sig,N_fft);

%% Blanchiment 
fft_y = medfilt1(abs(fft_sig), N_median);

%% Estimation
y = ifft(fft_y,N);
A = lpc(y,12);
out = filter(1,[0-A(2:end)],sig);
fft_out = fft(out.*win,N_fft);

figure
plot(f_vec, db(abs(fftshift(fft_out))))
hold on
plot(f_vec, db(abs(fftshift(fft_sig))))
xlim([0,fmax-df])

%% Différence entre les deux médianes glissantes
fft_y2 = medfilt1(abs(fft_out), N_median);

figure
plot(db(fft_y))
hold on
plot(db(fft_y2))
title('Médiane glissante avant et après blachimment')
xlabel('Fréquence')
ylabel('Amplitude (dB)')
legend('Original','Blanchit')
