function [p,h,f] = lotp(M,N,ws,alpha,taille_grid)
% function [p,h,f] = lotp(M,N,ws,alpha,taille_grid)
% cosine-modulated pour N pair
% p = filtre prototype
% h = filtres d'analyse
% f = filtres de synthese
% M = Nombre de bandes
% N = longueur des filtres
% ws = pi/(2*M) + epsilon
% alpha = Epb + alpha * Ebs
% taille_grid = nombre de points pour l'optimisation
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

tau = 0.5;
q=zeros(N/2,1);
err = 0.0001;
d=ones(taille_grid,1);
%Initialise p
inip=firls(N-1,[0. 1/(2*M)-1/(8*M) 1/(2*M)+1/(8*M) 1.],[1 1 0 0]);
p=inip(1:N/2);
wp=(0:pi/M/(taille_grid-1):pi/M)';
ind=(N-1)/2:-1:1/2;
C1=wp*ind;
CC1=2*cos(C1);
%Former Us
Us = zeros(N/2);
pimws =  (pi-ws)/2;
for i=1:N/2,
	arg1 = 2*i-N-1;
	for j=1:N/2,
		if(i == j),
			Us(i,i) = 4*(pimws-sin(arg1*ws)/(2*arg1));
		else
			imj = i-j;
			jmi=-imj;
			ipj = i+j-N-1;
			Us(i,j) = 4*(sin(imj*ws)/(2*jmi) - sin(ipj*ws)/(2*ipj));
		end
	end
end
dev = (p-q')*(p'-q);
iter = 0;
while(dev > err),
	iter = iter + 1
	Mp = CC1*p';
	Hp = diag(Mp);
	U = Hp*CC1;
	U = U + flipud(U);
	q=inv(U'*U+alpha*Us)*(U'*d);
	p=p-tau*(p-q');
	dev = sqrt((p-q')*(p'-q));
end;
p = [p fliplr(p)];
p=p';
n=((0:length(p)-1)-(N-1)/2)';
h = zeros(N,M);
for l=1:M,
	h(:,l)= 2*p.*cos((2*(l-1)+1)*pi*(n/(2*M)+0.25));
end;
f=conj(flipud(h));
