function [y,Vm,H]=CM_filterbank(h,g,M,x)
% Cosine-Modulated filter bank (Fig. 10.20)
% Input
%  h : FIR LPF coefficients or impulse response for analysis
%  g : FIR LPF coefficients or impulse response for synthesis
%  M : Number of channels (analysis-synthesis filter pairs)
%  x : Input sequence
%
% Output
%  y : Output sequence (delayed by M samples cod with the input x)
%  Vm: M-subband intermediate sequence matrix
%  H : MxN coefficient matrix of cosine-modulated filter bank
%
% Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only

if nargin < 4
    x = M;
    M = g;
    g = h;
end

%% Initialize
N = length(h);
D = N - 1;
m1 = [0 : M-1];
n1 = [0 : N-1];

X_buf = zeros(M,N);
W_buf = zeros(M,N);

Vm=[]; 

%% Creating filter bank
% Phase of modulating cosine
phis = repmat( (-1).^m1.'*pi/4 ,1,N);

% Instance of the modulating cosines for seperation and reconstruction
Cos  =2*cos(pi/M*(m1+0.5).'*(n1-D/2) + phis);
Cos_ =2*cos(pi/M*(m1+0.5).'*(n1-D/2) - phis);

% Initializing H (analysis) and G (synthesis) from the provided h et g
H = repmat(h,M,1).* Cos ;
G = repmat(g,M,1).* Cos_;

% X_buf and W_buf have the states of filters H and G.
lx = length(x); 
K  = ceil(lx/M); 
L  = K*M;

%% LOOOOOOOOOOOOOOOOOOOOOOOP
for n = 1:L+N-1
  if n<=lx,   
    xn=repmat(x(n),M,1);  
  else  
    xn=zeros(M,1);  
  end

  X_buf = [xn  X_buf(:,1:N-1)];
  Hx = sum(H.*X_buf,2); % Analysis filter bank outputs

  if mod(n-1,M)==0
    Vm = [Vm Hx]; 
    Wm = Hx; % Downsampling
  else
    Wm = zeros(M,1); % Upsampling
  end

  W_buf = [Wm  W_buf(:,1:N-1)];
  Gw = sum(G.*W_buf,2); % Synthesis filter bank outputs
  y(n) = sum(Gw);
end
