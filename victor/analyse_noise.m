N_FFT = 2^14;
FMAX = Fs/2;

N_WIN = 2*4096;
FLAG_A = 3200;
SIG = r(FLAG_A:FLAG_A+N_WIN-1);

FFT_SIG = fft(SIG.*hann(N_WIN), N_FFT);
FMAX = Fs/2;
dF = Fs/N_FFT;

F_VEC = [-FMAX + dF : dF : FMAX];
plot(F_VEC, db(abs(FFT_SIG)))
xlim([0 2000])
