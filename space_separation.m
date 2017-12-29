function [] = space_separation(input_filepath, out_dir, out_file)
% Perform subspace separation method on the file located at 'input_filepath',
% and writes the audio file at 'output_filepath'
%
% Arguments:
%     - input_filepath: file to perform separation on
%     - output_filepath: where to write the results to
%
% Returns:
%     -- None --

% Includes DESAM toolbox
startup

% Import and normalizing
[sig, Fs] = audioread(input_filepath);
sig = sig/max(1.1*abs(sig));
sig = extract_attack(sig, Fs);

% Subspace tracking method
[z, a, x] = analyse(sig, Fs);

% Separation
if length(x) < length(sig)
    r=x-sig(1:length(x));
else
    r=s-x(1:length(sig));
end

% Output filename
output_sig = strcat(out_dir,'sinus_',out_file);
output_noi = strcat(out_dir,'noise_',out_file);

% Data write
x = x/max(1.1*abs(x));
r = r/max(1.1*abs(r));

audiowrite(output_sig, x, Fs);
audiowrite(output_noi, r, Fs);

end
