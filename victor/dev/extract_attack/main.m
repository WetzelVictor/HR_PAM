[sig, Fs] = audioread('test.wav');
truncated_signal = extract_attack(sig, Fs);

figure
subplot(211)
plot(sig)
title('Original signal')

subplot(212)
plot(truncated_signal)
title('Truncated signal')
