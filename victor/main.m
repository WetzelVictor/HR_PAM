% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


[Micro, Fs] = audioread('data/Mesures_18_12_2017_matin/01-171218_plectredown_P1.wav');
[Acc, Fs] = audioread('data/Mesures_18_12_2017_matin/02-171218_plectredown_P1.wav');

deb = 1.351e5;
int = Fs/2;
s = Micro(deb:deb+int);


[z, a, x] =   analyse(s, Fs);

if length(x) < length(s)
    r=x-s(1:length(x));
else
    r=s-x(1:length(s));
end

audiowrite('output/sin_plectreP1.wav',x, Fs);
audiowrite('output/noise_plectreP1.wav',r, Fs);

%% 
N_FFT = 2^16;
N = length(r);
win = hann(N);
dF = Fs/N_FFT;
fmax = Fs/2;
f_vec = [-fmax+dF : dF:fmax];

fft_r = fft(r.*win,N_FFT);

figure
plot(f_vec, db(fftshift(abs(fft_r)))')
xlim([0,2000])
ylabel('Amp (dB)')
xlabel('Freq (Hz)')

