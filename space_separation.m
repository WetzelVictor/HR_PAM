function [] = space_separation(input_filepath, out_dir, out_file, fig_name)
% Perform subspace separation method on the file located at 'input_filepath',
% and writes the audio file at 'output_filepath'.
%
% Arguments:
%     - input_filepath: file to perform separation on
%     - output_filepath: where to write the results to
%
% Returns:
%     -- None --

% Includes DESAM toolbox
startup

% Import and normalizing
[sig, Fs] = audioread(input_filepath);
sig = sig/max(1.1*abs(sig));
sig = extract_attack(sig, Fs);

% Subspace tracking method
[z, alpha, x] = analyse(sig, Fs);
delta = log(abs(z));
f = 1/(2*pi)*angle(z);
a = abs(alpha);
phi = angle(alpha);
savedfile = replace(fig_name, '.fig', '.mat');
save(savedfile, 'delta', 'f', 'a', 'phi')
fmax = 15;
fen = 15;
fk = f(1:fmax, fen)*Fs;
ak = a(1:fmax, fen);
stem(fk,100+20*log(ak))
xlabel('Fréquence (Hz)')
ylabel('Amplitude(dB)')
savedfig = replace(fig_name, '.fig', '_params.fig');
savefig(savedfig)

% Separation
if length(x) < length(sig)
    r=x-sig(1:length(x));
else
    r=s-x(1:length(sig));
end

% Output filename
output_sig = strcat(out_dir,'sinus_',out_file);
output_noise = strcat(out_dir,'noise_',out_file);

% Data write
x = x/max(1.1*abs(x));
r = r/max(1.1*abs(r));

audiowrite(output_sig, x, Fs);
audiowrite(output_noise, r, Fs);

%% Spectrogram
load colors;

% Analysis parameters
N_fft = 2^15;
N_win = 512;
win = hann(N_win, 'periodic');
N_over = floor(N_win/8);


% AFFICHAGE DU SPECTROGRAMME
handles = figure;
spectrogram(r,win, N_over, N_fft,Fs,'yaxis','MinThreshold',-115)
colormap(colors)
xlabel('Temps (secondes)')
ylabel('Fr�quence (Hz)')
title(strcat('R�siduel ',out_file))
savefig(fig_name)
close all

end
