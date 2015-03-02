clear all;

rho = 0.9;
f = -0.5:0.01:0.5;
fh = 0:0.001:0.5;

% The different bit rates to be used
H1 = 0.75;
H2 = 2;
H3 = 5;

% Quantization noise
Deltasq1 = (2*pi*exp(1))/2^(2*H1);
Deltasq2 = (2*pi*exp(1))/2^(2*H2);
Deltasq3 = (2*pi*exp(1))/2^(2*H3);
sigmaq1 = Deltasq1/12;
sigmaq2 = Deltasq2/12;
sigmaq3 = Deltasq3/12;

% Input signal PSD
Sx = @(f) (1-rho^2)./(1+rho^2-(2*rho.*cos(2*pi.*f)));

% Input signal spectrum
Xf = sqrt(Sx(f));

%Kernel functions for the integrals that go into the Lagrange multiplier
IntSx1 = @(f) sqrt((sigmaq1.*Sx(f)));
IntSx2 = @(f) sqrt((sigmaq2.*Sx(f)));
IntSx3 = @(f) sqrt((sigmaq3.*Sx(f)));

%Perform integrals
q1 = integral(IntSx1,-0.5,0.5);
q2 = integral(IntSx2,-0.5,0.5);
q3 = integral(IntSx3,-0.5,0.5);

% Calculate Lagrange multipliers
Lagrange1 = (q1/(1+sigmaq1))^2;
Lagrange2 = (q2/(1+sigmaq2))^2;
Lagrange3 = (q3/(1+sigmaq3))^2;

% Matched receive filters
Hfsq1 = @(f) sqrt((Lagrange1.*Sx(f))/sigmaq1)-Lagrange1;
Hfsq2 = @(f) sqrt((Lagrange2.*Sx(f))/sigmaq2)-Lagrange2;
Hfsq3 = @(f) sqrt((Lagrange3.*Sx(f))/sigmaq3)-Lagrange3;

figure(1);
plot(fh,Hfsq1(fh),fh,Hfsq2(fh),fh,Hfsq3(fh));
title('Frequency response of matched receive filters');
xlabel('Normalized frequency (f)');
ylabel('H(f)');
legend('0.75 bits/sample','2 bits/sample','5 bits/sample');

% Matched transmit filters
Gfsq1 = @(f) sqrt(sigmaq1./(Lagrange1.*Sx(f)))-(sigmaq1./Sx(f));
Gfsq2 = @(f) sqrt(sigmaq2./(Lagrange2.*Sx(f)))-(sigmaq2./Sx(f));
Gfsq3 = @(f) sqrt(sigmaq3./(Lagrange3.*Sx(f)))-(sigmaq3./Sx(f));

figure(2);
plot(fh,Gfsq1(fh),fh,Gfsq2(fh),fh,Gfsq3(fh));
title('Frequency response of matched transmit filters');
xlabel('Normalized frequency (f)');
ylabel('G(f)');
legend('0.75 bits/sample','2 bits/sample','5 bits/sample');

% Output noise PSD
Nf1 = sigmaq1.*Hfsq1(fh);
Nf2 = sigmaq2.*Hfsq2(fh);
Nf3 = sigmaq3.*Hfsq3(fh);

% Output signal PSD
Yf1 = Sx(fh).*Hfsq1(fh).*Gfsq1(fh);
Yf2 = Sx(fh).*Hfsq2(fh).*Gfsq2(fh);
Yf3 = Sx(fh).*Hfsq3(fh).*Gfsq3(fh);

figure(3);
semilogy(fh,Yf1,'c');
hold on;
semilogy(fh,Nf1);
semilogy(fh,Yf2,'m');
semilogy(fh,Nf2,'g');
semilogy(fh,Yf3,'k');
semilogy(fh,Nf3,'r');
title('Output signal and noise power spectrums');
xlabel('Normalized frequency (f)');
ylabel('Y(f), N(f) (logarithmic scale)');
legend('Signal (0.75 bits/sample)','Noise (0.75 bits/sample)','Signal (2 bits/sample)','Noise (2 bits/sample)','Signal (5 bits/sample)','Noise (5 bits/sample)');
hold off;

% Signal-to-Noise ratios
SNR1 = 10*log10(sum(Yf1)/sum(Nf1));
SNR2 = 10*log10(sum(Yf2)/sum(Nf2));
SNR3 = 10*log10(sum(Yf3)/sum(Nf3));