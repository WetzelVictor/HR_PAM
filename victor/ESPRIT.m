function [dk, fk] = ESPRIT(x, n, K)
% Méthode HR ESPRIT: graphique
% INPUT:
%   x: signal
%   n: longueur vecteur de données
%   K: dimension de l'espace signal

%% INIT
N = length(x);
l = N - n + 1;

% Matrice X (matrice de Hankel)
X = hankel(x(1:n), x(n:N));

% Matrice de corrélation
RxxB_mat = 1/l * X * X';

%% ESTIMATION ESPACE SIGNAL
[U1, vp_Rxx, U2] = svd(RxxB_mat);

W = U1(:,1:K);

%% ESTIMATION FREQ & F.AMORTI.
Wdo = W(1:n-1,:);
Wup = W(2:n,:);

% Estimation
Phi = pinv(Wdo) *Wup;
vp_Phi = eig(Phi); % estimation des pôles complexes

% delta k et fk
dk = log(abs(vp_Phi)); 
fk = (1/(2*pi)) * angle(vp_Phi);
end
