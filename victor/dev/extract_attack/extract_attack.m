function [truncated_sig] = extract_attack(sig, Fs)
% extract the transient from a clean signal (e.g. where the attack emerges from
% the rest of the sound scene. The method is a simple search-for-the-maximum
% type
%
% Arguments:
%   - signal: the vector you want to extract the attack from
%   - Fs: frequency samplerate
%
% Returns:
%   - truncated_sig : truncated signal

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
