function [dk, fk] = MUSIC(x, n, K)
% Méthode HR MUSIC: graphique
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

W = U1(:,K+1:end);

%% Méthode MUSIC
f = [0:0.001:1];
d = [-0.1:0.001:0.1]';

for k1=1:length(f)
  for k2=1:length(d)
    zk = exp(d(k2) + 1j*2*pi*f(k1));

    lognu = [0:n-1]*log(zk);
    nu = exp(lognu).';

    P(k1,k2) = 1/(norm(W'*nu).^2);
  end
end

mesh(20*log(P));
zlabel('Amplitude')
ylabel('Fréquence (Hz)')
xlabel('\delta')

end
