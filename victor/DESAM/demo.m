% function demo()
%
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


[s, Fs] = audioread('ClocheB.wav');
s = s/max(abs(s));
% s = s(1:85000);
s = s + randn(size(s))*0.000001;

[z, a, x] = analyse(s, Fs);

if length(x) < length(s)
    r=x-s(1:length(x));
else
    r=s-x(1:length(s));
end

audiowrite('highResolution_telecomParisTech_sin.wav',x, Fs);
audiowrite('highResolution_telecomParisTech_noise.wav',r, Fs);
specgram(r, 1024,Fs)
title('Résiduel')
figure
specgram(x, 1024,Fs)
title('Sinusoïdes')
