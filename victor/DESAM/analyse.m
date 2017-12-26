function [The_z,The_alpha,x] = analyse(s,Fs)
%[z,alpha,x] = analyse(s,Fs) effectue l'analyse du signal x
%z contient les p?les complexes estim?s
%alpha contient les amplitudes complexes estim?es
%x est le signal resynth?tis?
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

% DEFINITION DU BANC DE FILTRES D'ANALYSE
D = 3; % facteur de d?cimation dans chaque sous-bande = nombre de sous-bandes de fr?quences positives
M = 2*D; % nombre total de sous-bandes (D sous-bandes de fr?quences positives + D sous-bandes de fr?quences n?gatives)
L = 128; % longueur des filtres d'analyse
h = filterbank(L,D); % calcul des filtres d'analyse


% INITIALISATIONS
s2 = filter([1,-0.98],1,s); % pr?-accentuation
n = 64; % dimension des vecteurs de donn?es
l = 65; % nombre de vecteurs par fen?tre d'analyse
N = n+l-1; % longueur totale de l'horizon d'observation
The_Freq = []; % Matrice contenant les fr?quences estim?es ? chaque instant d'analyse
The_z = []; % Matrice contenant les p?les complexes estim?s ? chaque instant d'analyse
The_alpha = []; % Matrice contenant les amplitudes complexes estim?es ? chaque instant d'analyse
rang = [20,20,20]; % nombre de fr?quences recherch?es dans les sous-bandes ; d?tection automatique si 0


% ANALYSE DU SIGNAL
for k=1:length(rang), % pour toutes les sous-bandes
   y = filter(h(:,k),1,s2); % filtrage
   y = y(1:D:end); % d?cimation
   [z,alpha] = estim(y,rang(k),n,l); % analyse HR adaptative
   Freq = mod(angle(z)/(2*pi),1); % calcul des fr?quences dans les sous-bandes
   % INITIALISATIONS
   Freq2 = zeros(size(Freq)); % fr?quences en bande pleine
   z2 = zeros(size(z)); % p?les complexes en bande pleine
   alpha2 = zeros(size(alpha)); % ampltiudes complexes en bande pleine
   % REMAPPAGE DES FREQUENCES EN BANDE PLEINE
   for t=1:size(Freq,2), % pour tous les instants d'analyse
		if mod(k,2) ==1, % cas des sous-bandes impaires
		   I = Freq(:,t)<.5; % on ne conserve que les fr?quences comprises entre 0 et 0.5 (les autres sont alias?es)
			Freq2(1:sum(I),t) = ((k-1)+2*Freq(I,t))/M; % remappage
		else % cas des sous-bandes paires
		   I = Freq(:,t)>=.5;  % on ne conserve que les fr?quences comprises entre 0.5 et 1 (les autres sont alias?es)
			Freq2(1:sum(I),t) = ((k-1)+2*(Freq(I,t)-.5))/M; % remappage
      end;
      z2(1:sum(I),t) = exp(log(abs(z(I,t)))/D+i*2*pi*Freq2(1:sum(I),t)); % remappage des p?les en bande pleine
      % CORRECTION DES AMPLITUDES COMPLEXES
      V = exp(log(1./z2(1:sum(I),t))*(0:size(h,1)-1)); % matrice de Vandermonde construite sur les p?les estim?s
      alpha2(1:sum(I),t) = alpha(I,t) ./ ( (V*h(:,k)) .* (V(:,1:2)*[1;-0.98]) ) ; % compensation du filtre de pr?-accentuation
   end;
	% JUXTAPOSITION DES PARAMETRES ESTIMES A PARTIR DES SOUS-BANDES  
   The_Freq = [The_Freq;Freq2];
   The_z = [The_z;z2];
	The_alpha = [The_alpha;alpha2];
end;

% RESYNTHESE DU SIGNAL
x = 2*real(synthese(The_z,The_alpha,D,n,l,1));

% % AFFICHAGE DU HR-OGRAMME
load colors; % chargement des couleurs pour le HR-ogramme
% The_Pow = abs(The_alpha).^2/n; % calcul des puissances estim?es
% seuil = 150; % seuil en dB pour r?gler le contraste du HR-ogramme
% HRogram(N*D/8,Fs,The_Freq,10*log10(The_Pow),1,seuil,colors); % affichage

% AFFICHAGE DU SPECTROGRAMME
figure; specgram(s,1024,Fs);
colormap(colors);
xlabel('Temps (secondes)');
ylabel('Fr?quence (Hz)');
