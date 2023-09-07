function [signal] = sigInj(fpos, ta, phase, fmin, fmax, m1,m2,r,datalen, initial_phase, N, avec, A)
%SIGINJ Returns time domain vector of waveform for injection with custom
%strain amplitude
%   This function creates a time domain vector of a custom injected
%   waveform normalized to a strain amplitude which is a function of the
%   component masses and the distance to the source. 
% This factor is given in Cutler and Flanagan 1994 for the newtonian waveform 
% (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.49.2658)
%  Input: fpos: Positive frequency vector, 
%         ta: time of arrival, 
%         phase: coalescence phase, 
%         [fmin, fmax]: waveform frequency bounds, 
%         {m1,m2,r}: the component masses and the distance in Mpc respectively, 
%         datalen: data length in seconds, 
%         initial_phase: initial phase, 
%         N: total number of samples, 
%         avec: Precalculated vectors for phase generation, 
%         A: frequency magnitude vector
% Output: signal: time-domain vector of waveform with normalized strain amplitude 

%Constants
% fmin = 30;
c = 3*10^8;
G = 6.6743*10^-11;
Msolar = 1.989*10^30;
Mpc = 3.8057*10^22; %1 Megaparsec in meters
% m1_val = m1;
% m2_val = m2;
m1_val = m1*Msolar;
m2_val = m2*Msolar;
M = m1_val + m2_val;
u = m1_val*m2_val/(m1_val + m2_val);
chirpmass = (u^3*M^2)^(1/5);
% pzi = 3*((chirpmass/Msolar)^(-5/3))*(fmin/100)^(-8/3);
% Nfac = (1/sqrt(50))*(fmin/100)^(2/3);
r = r*Mpc;
%Factors to keep phase vector unnormalized
snr = 1;
normfac = 1;

%Initiate parameters
% zhi = 34.54*((chirpmass/Msolar)^(-5/3));
% 
% Afac = (1.92*10^(-23))*((zhi/25)^(-1))*((r/(100*Mpc))^(-1));
% 
% Amplfac = Afac*sqrt(zhi)*((2/(3*fmin))^(1/2))*((1/fmin)^(-7/6));

%Multiply strain amplitude factor

% Afac = ((2*u*G)/(r*c^4))*(G*M*pi)^(2/3);
Afac = ((384/5)^(1/2))*(pi^(2/3))*(chirpmass^(5/6))*(G^(5/6))*(c^(-3/2))/r;

%Create Fourier Phase vector
wavephase = gen2PNwaveform(fpos, ta, phase, fmin, fmax, m1,m2,datalen, initial_phase, snr, N, avec, normfac);
% wavephase = gen2PNwaveform_tau(fpos, ta, phase, fmin, fmax, m1,m2,datalen, initial_phase, snr, N, avec, normfac);

%Multiply Fourier Magnitude vector and strain amplitude factor
wavefourier = Afac*((1/fmin)^(-7/6))*A.*wavephase;
% wavefourier = ((1/fmin)^(-7/6))*A.*wavephase;

%Create waveform vector in time domain
signal = ifft(wavefourier);

end

