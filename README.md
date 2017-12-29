# HR_PAM
PAM project, part of the ATIAM Master at Jussieu-IRCAM-Télécom ParisTech

## Signal processing
Subspace base method for instrument making property extraction.

## Systematic extraction: Manual
You have to install the matlab engine module, with a default 64bit python install
(not anaconda install). Instructions can be found on this
[link](https://fr.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html).


First, you need to download the audio files from the project, and move them to
the `data` folder. Then run the `replicate_datastructure.py` script. Finally,
execute `sys_separation.py` to perform noise/signal separation.
