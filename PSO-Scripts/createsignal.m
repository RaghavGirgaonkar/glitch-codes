function [signal] = createsignal(siglen, frange, sampFreq, masses, r, initial_phase, phase, ta, mfac)
%Signal Injection code

%Create Positive DFT frequency vector
N = siglen*sampFreq;
fpos = (0:floor(N/2))*(1/siglen);

%Get frequency range
fmin = frange(1);
fmax = frange(2);

%Preprocessing
[A,avec, ~] = preprocessing(fmin,fmax,fpos, siglen, N);

%Create Signal
m1 = masses(1);
m2 = masses(2);

signal = mfac*sigInj(fpos, ta, phase, fmin, fmax, m1,m2,r,siglen, initial_phase, N, avec, A);
end

