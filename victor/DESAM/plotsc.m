function plotsc(x,y,z,colors);
%plotsc(x,y,z,colors) affiche le point de coordonnées (x,y) dans le plan avec une couleur de niveau z dans la palette colors
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

N = length(x);
NCol = size(colormap,1);
dx = x(2)-x(1);
if (length(y) ~= N) | (length(z) ~= N),
   disp('??? Error using ==> plotsc');
   disp('Vectors must be the same lengths.');
elseif max(abs(imag(z))) ~=0,
   disp('??? Error using ==> plotsc');
   disp('Third argument must be real.');
else
   hold on;
   for k=1:N,
		h = plot([x(k),x(k)+dx],[y(k),y(k)]);
      c = z(k);
      set(h,'linewidth',[1]); % 20
	   set(h,'color', colors(ceil(c*NCol+eps),:)); % [c,c,c]);
   end;
   hold off;
end;