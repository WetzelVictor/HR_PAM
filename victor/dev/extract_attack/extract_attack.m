function [truncated_sig] = extract_attack(sig, Fs)

% Retrieving maximum and its index
[threshold, index] = max(abs(sig));

% Length before 
len = 0.8; % length (sec.)
len_sample = len*Fs;

% Retrieving truncated signal
flagA = index - len_sample;
flagB = index + 4*len_sample;

truncated_sig = sig(flagA: flagB);

end
