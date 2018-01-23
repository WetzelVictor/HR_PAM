function [] = plot_params(inputpath, filename)

load(inputpath);
Fs = 48000;
fmax = 40;
fen = 15;
fk = f(1:fmax, fen)*Fs;
ak = a(1:fmax, fen)*Fs;

stem(fk, 20*log(ak))
xlabel('Fréquence (Hz)')
ylabel('Amplitude(dB)')
savefig(filename)

end
