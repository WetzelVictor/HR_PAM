function [] = plot_modes_struct(inputpath, filename, fig_name)

load(inputpath);
Fs = 48000;
tmax = 50;
n = 4;
p = size(f,1);
data = zeros(p,2,tmax);
modA = zeros(1,tmax);
modT = zeros(1,tmax);
modF = zeros(n,1,tmax);

for k = 1:tmax
    data(:,:,k) = [f(:,k)*Fs, a(:,k)*Fs];
    data(:,:,k) = sortrows(data(:,:,k),1);
    fk = data(:,1,k);
    ind = find(fk,n);
    for i = 1:n
        if fk(ind(i)) <= 180 && fk(ind(i)) >= 80
            modA(k) = data(ind(i),2,k);
        elseif fk(ind(i)) <= 300 && fk(ind(i)) >= 210
            modT(k) = data(ind(i),2,k);
        end
    end 
    modF(:,k) = data(ind,1,k);
end

%  disp(modF)
stem(20*log(modA), 'b')
hold on
stem(20*log(modT), 'r')
hold off
xlabel('Temps (frames)')
ylabel('Amplitude(dB)')
legend('A0','T1')
title(filename)
savefig(fig_name)

end