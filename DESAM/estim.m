function [z,alpha] = estim(y,r,n,l);
% [z,alpha] = estim(y,r,n,l) effectue l'analyse HR adaptative du signal y
% r est le nombre de fréquences recherchées
% si r=0, ce nombre est détecté automatiquement
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

% INITIALISATIONS
N = n+l-1; % longueur de l'horizon d'observation
N8 = floor(N/8);
N4 = 2*N8;
rmax = n/2; % valeur maximale admissible du paramètre r
q = 51; % longueur du filtre de rang pour éliminer les pics dans le périodogramme
p = 5; % longueur du filtre de blanchiment du bruit
y = [zeros(p-1,1);y(:)]; % insertion de zeros au début du signal
duree = p+l-1:N8:length(y)-n+1; % instants d'analyse
a = zeros(p,length(duree)); % filtre de blanchiment
win_N = hanning(N); % utilisation d'une fenêtre de Hanning pour le calcul du périodogramme
e = zeros(p,1);
e(1) = 1;
Np = 2^nextpow2(2*N-1);
ax = (0:Np-1)/Np;
k = 1; 

% BLANCHIMENT DU BRUIT ET SELECTION DE L'ORDRE
for t = duree; % pour chaque fenêtre d'analyse
   x = y(t-l+1:t+n-1); % extraction d'un segment de signal
   P = abs(fft(x.*win_N, Np)).^2/N; % calcul du périodogramme
   % LISSAGE DU PERIODOGRAMME PAR FILTRE DE RANG
   H = hankel([P(end-(q-1)/2+1:end);P(1:(q+1)/2)],[P((q+1)/2:end);P(1:(q-1)/2)]);
   H = sort(H,1);
   P = H(q/3,:).';
   % PREDICTION LINEAIRE DU BRUIT
   rx = ifft(P);
   a(:,k) = toeplitz(rx(1:p),conj(rx(1:p))) \ e;
   sigma2 = abs(a(1,k));
   a(:,k) = a(:,k) / sigma2;
   if r==0, % s'il faut estimer l'ordre du modèle
	   % BLANCHIMENT DU BRUIT
	   x = filter(a(:,k),1,y(t-l-p+2:t+n-1)); % filtrage
	   x = x(p:end); % élimination du transitoire
	   % CALCUL DES PRINCIPAUX VECTEURS PROPRES
	   H = hankel(x(1:n),x(n:N));
	   if k==1, % à l'initialisation, calcul exact
			[U,S,V] = svd(H,0);
	      U = U(:,1:rmax);
	   else, % ensuite, algorithme adaptatif
	      A = H*(H'*U);
	      [U,R] = qr(A,0);
	   end;
	   % CALCUL DU CRITERE ESTER
	   W = zeros(n,0);
		Xi = zeros(n-1,0);
	   Psi = zeros(0,0);
	   for r2 = 1:rmax,
			wr = U(:,r2);
			wl = U(end,1:r2)';
			psir = W(1:end-1,:)'*wr(2:end);
			psil = W(2:end,:)'*wr(1:end-1);
			psilr = wr(1:end-1)'*wr(2:end);
			Psi = [Psi,psir;psil',psilr];
			phi = Psi' * wl;
			xi = wr(2:end) - W(1:end-1,:)*psir - wr(1:end-1)*psilr;
			W = U(:,1:r2);
	      Xi = [Xi - wr(1:end-1)*psil',xi];
	      if abs(1 - wl'*wl)>eps,
	         E = Xi - (W(1:end-1,:)*wl*phi') / (1 - wl'*wl);
	         ESTER(k,r2) = norm(E)^2;
	      end;
	   end;
	   % SELECTION DE L'ORDRE r
	   if sum(ESTER(k,:))>0,
	      rang(k) = max((ESTER(k,:)<0.01).*(1:rmax));
	   else
	      rang(k) = 0;
      end;
   end; % fin de sélection de l'ordre
   k = k+1; % itération suivante
end;

if r==0, % s'il faut estimer l'ordre du modèle
	% FUSION DES INFORMATIONS OBTENUES A CHAQUE INSTANT D'ANALYSE ET SELECTION D'UN ORDRE r GLOBAL
	rang = sort(rang);
	r = rang(round(length(rang)*9/10));
end;

% ANALYSE HR ADAPTATIVE
k = 1;
t1 = N4+1;
t2 = N/4+1;
for t = duree; % pour chaque fenêtre d'analyse
   % BLANCHIMENT DU BRUIT
   x = filter(a(:,k),1,y(t-l-p+2:t+n-1));
   x = x(p:end);
   % CALCUL DES PRINCIPAUX VECTEURS PROPRES
   H = hankel(x(1:n),x(n:N));
   if k<=2, % à l'initialisation, calcul exact
      [U,S,V] = svd(H);
      W = U(:,1:r);
   else, % ensuite, algorithme adaptatif
      A = H*(H'*W);
      [W,R] = qr(A,0);
   end;
   % METHODE ESPRIT
   Phi = pinv(W(1:end-1,:)) * W(2:end,:);
   z(:,k) = eig(Phi);
   % CALCUL DE LA MATRICE DE VANDERMONDE NORMALISEE
   V = zeros(n+l-1,r);
   for k2=1:r,
      if abs(z(k2,k)) < 1,
   		V(:,k2) = z(k2,k).^((0:n+l-2)');
	   else,
	      V(:,k2) = z(k2,k).^((-(n+l-2):0)');
   	end;
	   V(:,k2) = V(:,k2) / norm(V(:,k2));
   end;

% ESTIMATION DES AMPLITUDES COMPLEXES
alpha(:,k) = pinv(V)*x; % moindres carrés
V = exp(log(1./z(:,k))*(0:p-1)); % matrice de Vandermonde
alpha(:,k) = alpha(:,k) ./ (V*a(:,k)); % correction pour tenir compte du filtre de blanchiment
k = k+1; % itération suivante
t1 = t1 + N8; % recouvrement de 7/8 entre deux fenêtres consécutives
t2 = t2 + N/8;

end;
