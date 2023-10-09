function [signal] = createsignal(siglen, frange, sampFreq, masses, r, initial_phase, ta, mfac, theta,phi, psi)
%Signal Injection code
%mfac = 3.2 for GW170817 SNR of 26 in Livingston file 

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

hplus = mfac*sigInj(fpos, ta, 0, fmin, fmax, m1,m2,r,siglen, initial_phase, N, avec, A);
hcross = mfac*sigInj(fpos, ta, pi/2, fmin, fmax, m1,m2,r,siglen, initial_phase, N, avec, A);

Fplus = ((1 + cos(theta)^2)/2)*cos(2*phi)*cos(2*psi) - cos(theta)*sin(2*phi)*sin(2*psi);
Fcross = ((1 + cos(theta)^2)/2)*cos(2*phi)*sin(2*psi) +  cos(theta)*sin(2*phi)*cos(2*psi);

signal = Fplus*hplus + Fcross*hcross;
end

