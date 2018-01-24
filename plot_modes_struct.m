function [] = plot_modes_struct(inputpath, filename, fig_name)

load(inputpath);
Fs = 48000;
L = 16; % longueur trame en échantillons
tmax = 70; % trame temporelle max
n = 4; % nb de premiers modes considérés
p = size(f,1);
data = zeros(p,3,tmax);
modA = zeros(1,tmax); % tableau des modes acoustiques
modT = zeros(1,tmax); % tableau des modes de table
modF = zeros(n,1,tmax); % tableau des fréquences recherchées

for k = 1:tmax % pour chaque trame
    data(:,:,k) = [f(:,k)*Fs, a(:,k), delta(:,k)];
    data(:,:,k) = sortrows(data(:,:,k),1); % tri par f croissant
    fk = data(:,1,k);
    ind = find(fk,n); % recherche des n premiers f non nulles
    for i = 1:n
        fi = fk(ind(i));
        ai = data(ind(i),2,k);
        deltai = data(ind(i),3,k);
        if fi <= 168 && fi >= 162 % détection modes A pour plexiglas
            modA(k) = ai*exp(deltai*k*L);
%         elseif fi <= 280 && fi >= 270 % détection modes T pour bois
%             modT(k) = ai*exp(deltai*k*L);
        end
    end 
    modF(:,k) = data(ind,1,k);
end

% disp(modF)
% for i =1:n
%     plot(modF(i,:))
%     hold on
% end
% hold off
x = (1:tmax)*L/Fs*1000;
xa = linspace(0, 30);
ya = 10*log10(modA.^2);
yt = 10*log10(modT.^2);
inda = isfinite(ya);
indt = isfinite(yt);
pa = polyfit(x(inda), ya(inda), 1);
fa = polyval(pa, xa);
pt = polyfit(x(indt), yt(indt), 1);
ft = polyval(pt, xa);
scatter(x, ya, 'b', 'filled')
hold on
% scatter(x, yt, 'r', 'filled')
plot(xa, fa, 'b--')
% plot(xa, ft, 'r--')
hold off
axis ([0 tmax*L/Fs*1000 -80 0])
xlabel('Temps (ms)')
ylabel('Amplitude (dB)')
legend('m1 165 Hz')
title(filename)
savefig(fig_name)

end