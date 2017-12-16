% TP: méthode HR
% Victor Wetzel
% ATIAM 2017

clear all; close all;

%% INIT
% [sig, Fs] = audioread('clocheA.wav');

%% Exponentielles synthétiques
f0 = 1/4;     % Freq fondamentale f0
N = 63;
f1 = f0 + 1/N;

a0 = 1; a1 = 10;
a = [a0; a1];
d0 = 0; d1 = -0.05;
d = [d0; d1];
phi0 = rand(1,1)*2*pi; phi1 = rand(1,1)*2*pi;
phi = [phi0; phi1];

% Vecteur fréquence
f = [-1/2:1/N:1/2 - 1/N];

x = Synthesis(N, [d0 d1],[f0 f1],[a0 a1], [phi0 phi1]);

%% Méthode ESPRIT
n = 32; K = 2;
ESPRIT(x, n, K);

[a_est, phi_est] = LeastSquares(x,[d0;d1],[f0;f1]);

%% Méthode MUSIC
K = 20;
n = 32;
MUSIC(x,n,K);

%%
[x, Fe] = audioread('ClocheA.WAV');

[dk, fk] = ESPRIT(x,512,54);
[ak, phik] = LeastSquares(x, dk, fk);
y = Synthesis(100000,dk,fk,ak,phik);

