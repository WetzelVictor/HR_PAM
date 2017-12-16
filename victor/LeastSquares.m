function [a, phi] = LeastSquares(x, delta, fk)
% Méthode des moindes carrés pour estimer les amplitudes et phases
%
% INPUT:
%   x: signal
%   delta (vec): amortissements
%   f (vec): fréquences
%
% OUTPUT:
%   a: amplitudes
%   phiK: phases

%% INIT
N = length(x);
t = [1:N]';
k = length(fk);
vecPhiF = delta.' + 1j*2*pi*fk.';

% Création matrice de Van der Monde
preV = t*vecPhiF;
Vn = exp(preV);

alpha =  pinv(Vn)*x;

%% OUTPUT
a = abs(alpha);
phi= angle(alpha);

end
