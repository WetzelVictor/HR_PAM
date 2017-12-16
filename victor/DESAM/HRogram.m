function HRogram(N8,Fs,Freq,Pow,M,seuil,colors,Name);
% HRogram(N8,Fs,Freq,Pow,M,seuil,colors,Name) affiche le HR-ogramme du signal analysé
%
% Copyright (C) 2004-2008 Roland Badeau
% Send bugs and requests to roland.badeau@telecom-paristech.fr
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

% Si le titre de la figure n'est pas donné, on lui assigne une chaîne vide
if nargin < 8,
   Name = '';
end

h = get(gca);
hold on;
Freq = Freq(:,1:M:end); % décimation de la matrice des fréquences pour gagner du temps
Pow = Pow(:,1:M:end); % décimation de la matrice des puissances pour gagner du temps
Pow = max(Pow+seuil,0); % réglage du contraste
Pow = Pow ./ max(max(Pow)+eps); % normalisation des puissances
t = (0:size(Freq,2)-1)*M*N8/Fs; % instants d'analyse
h = rectangle('Position',[0,0,max(t),Fs*max(max(Freq))*1.1]); % affichage du fond
set(h,'FaceColor',colors(1,:)); % couleur du fond
for n=1:size(Freq,1),
   plotsc(t, Fs*Freq(n,:), Pow(n,:),colors);
end
hold off;
axis([0,max(t),0,Fs*max(max(Freq))*1.1]); % réglage des axes
title(Name);
xlabel('Temps (secondes)');
ylabel('Fréquence (Hz)');

