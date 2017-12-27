function h = filterbank(N,D);
%h = filterbank(N,D) renvoie les D filtres d'analyse de longueur N
%
% Copyright (C) Jacques Prado
% Send bugs and requests to jacques.prado@telecom-paristech.fr
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

% CALCUL DU FILTRE PROTOTYPE D'UN BANC DE FILTRES A DCT
lfft = 4096;
fg=0:1/lfft:1/2;
ws=0.85*pi/D+pi/1000;
alpha =100;
taille_grid=200;            % taille_grid = nombre de points pour l'optimisation
if (rem(N,2)==0),
	[p, h, f] = lotp(D,N,ws,alpha,taille_grid);
else
	[p, h, f] = loti(D,N,ws,alpha,taille_grid);
end
% TRANSLATIONS DU FILTRE PROTOTYPE POUR OBTENIR LES FILTRES D'ANALYSE DU BANC A TFD
h = zeros(N,D);
n=((0:length(p)-1)-(N-1)/2)';
p0= 2*p.*cos(pi*(n/(2*D)+0.25));
M = 2*D;
for l=1:D,
   h(:,l)= p0.* exp(2*pi*i*(l-.5)*n/M);
end;
