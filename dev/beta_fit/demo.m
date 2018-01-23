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
clear all
close all

load('exemple_separation.mat')
ind = 5;
z = z(:,ind);
a = a(:,ind);

%% Coefficients de qualité
f0 = Fs*mod(angle(z)/(2*pi),1);
dk = log(abs(z)); 
tau = 1./(2*Fs*dk);
Q = -2*pi*f0.*tau;
n = 1:200;
f00 = sort(f0(f0>20)); f00 = f00(1);
nF = n*f00;


flim = 3000;
Qs = 10000;
Qmin = -100;
Qmax = 1e8;

Nfreq = length(f0);

% INIT
% Freq
f_string = [];
f_body = [];
% Poles
z_string = [];
z_body = [];
% Qfactor
Q_string = [];
Q_body = [];
% index
ind_string = [];
ind_body = [];
% amplitude
a_body = [];
a_string = [];

distance = [];
for i=1:Nfreq
    % First condition
    minCombDistance = min(abs(nF - f0(i)));
    distance = [distance minCombDistance];
    
    if f0(i) < 20 
        continue   
    elseif f0(i) > flim
        f_string = [f_string f0(i)];
        Q_string = [Q_string Q(i)];
        ind_string = [ind_string i];
        z_string = [z_string z(i)];
        a_string = [a_string a(i)];
%     elseif Q(i) < Qmax % &&Q(i) > Qmin
    elseif minCombDistance < 10
            f_string = [f_string f0(i)];
            Q_string = [Q_string Q(i)];
            ind_string = [ind_string i];
            z_string = [z_string z(i)];
            a_string = [a_string a(i)];
        elseif Q(i) < Qs
            f_body = [f_body f0(i)];
            z_body = [z_body z(i)];
            Q_body = [Q_body Q(i)];
            ind_body = [ind_body i];
            a_body = [a_body a(i)];

%         end
    end
end

figure
stem(f_string,ones(size(f_string)))
xlim([0 3000])
ylim([0.8 1.2])


Nf = length(f_body);
t = [0:1:100000];
imp = zeros(size(t));
for i = 1:Nf
    imp = imp + a_body(i)*z_body(i).^t;
end

imp = real(imp(:));
Nfft = length(imp)*2;
win = hann(length(imp));
IM = fft(imp.*win, Nfft);
f= [0:1/Nfft:1-1/Nfft]*Fs;
figure
plot(f,db(abs(IM)))
xlim([20 3000])
    
    