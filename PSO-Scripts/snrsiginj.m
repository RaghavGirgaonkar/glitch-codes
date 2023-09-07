function [signal] = snrsiginj(psd, snr, masses, frange, ta, phase, datalen, sampFreq, initial_phase)
%Inject signal of user specified SNR

%% Constants
c = 3*10^8;
Msolar = 1.989*10^30;
G = 6.6743*10^-11;

%Create Positive DFT frequency vector
N = datalen*sampFreq;
fpos = (0:floor(N/2))*(1/datalen);

%Get frequency range
fmin = frange(1);
fmax = frange(2);

%Preprocessing
[A,avec, ~] = preprocessing(fmin,fmax,fpos, datalen, N);

%Create Signal
m1 = masses(1);
m2 = masses(2);

%Create Normalization factor
TF = 1./sqrt(psd);
%Set Values of TF below fmin to 0
TF(1:fmin*datalen) = 0;

%Create entire Transfer Function vector
negFStrt = 1-mod(N,2);
kNyq = floor(N/2)+1;

TFtotal = [TF, TF((kNyq-negFStrt):-1:2)];

%Create AbysqrtPSD
AbysqrtPSD = A.*TFtotal;

% Scalar factor of 1/N is due to Parseval's theorem
dataLen = N;
innProd = (1/dataLen)*(AbysqrtPSD)*AbysqrtPSD';
genNormfacSqr = real(innProd);
genNormfac = 1/sqrt(genNormfacSqr);

%Create nomralized signal with specified SNR
m1_val = m1*Msolar;
m2_val = m2*Msolar;
M = m1_val + m2_val;
u = m1_val*m2_val/(m1_val + m2_val);
n = u/M;
tau0 = (5/(256*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-5/3))*(1/n);
tau1p5 = (1/8)*(1/fmin)*((G*M*pi*fmin/c^3)^(-2/3))*(1/n);

% wavephase = gen2PNwaveform_tau(fpos, ta, phase, fmin, fmax, tau0,tau1p5,datalen, initial_phase, snr, N, avec, genNormfac);
wavephase = gen2PNwaveform(fpos, ta, phase, fmin, fmax, m1,m2,datalen, initial_phase, snr, N, avec, genNormfac);
wavefourier = A.*wavephase;

%Create waveform vector in time domain
signal = ifft(wavefourier);
end

