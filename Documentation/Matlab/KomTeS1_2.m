clear all;

NL = 29; % Number of points for low resolution spectrum
NH = 999; % Number of points for high resolution spectrum
fbh = -0.5:1/NH:0.5; % bilateral half spectrum frequency vector
fuh = 0:1/NH:0.5; % unilateral half spectrum frequency vector
flf = 0:1/NL:1-1/NL; % Low resolution unilateral full spectrum frequency vector
fli = 1:ceil(length(flf)/2); % Half spectrum index vector
fhf = 0:1/NH:1-1/NH; % High resolution unilateral full spectrum frequency vector
fhi = 1:ceil(length(fhf)/2); % Half spectrum index vector

% The different bit rates to be used
H1 = 0.75;
H2 = 2;
H3 = 5;
rho = 0.9;

%------------ Exercise 2a ------------%

% Quantization noise
Deltasq1 = (2*pi*exp(1))/2^(2*H1);
Deltasq2 = (2*pi*exp(1))/2^(2*H2);
Deltasq3 = (2*pi*exp(1))/2^(2*H3);
sigmaq1 = Deltasq1/12;
sigmaq2 = Deltasq2/12;
sigmaq3 = Deltasq3/12;

% Input signal PSD
Sx = @(f) (1-rho^2)./(1+rho^2-(2*rho.*cos(2*pi.*f)));

%Kernel functions for the integrals that go into the Lagrange multiplier
IntSx1 = @(f) sqrt((sigmaq1.*Sx(f)));
IntSx2 = @(f) sqrt((sigmaq2.*Sx(f)));
IntSx3 = @(f) sqrt((sigmaq3.*Sx(f)));

%Perform integrals
%q1 = integral(IntSx1,-0.5,0.5);
%q2 = integral(IntSx2,-0.5,0.5);
%q3 = integral(IntSx3,-0.5,0.5);

% Calculate Lagrange multipliers
Lagrange1 = (integral(IntSx1,-0.5,0.5)/(1+sigmaq1))^2;
Lagrange2 = (integral(IntSx2,-0.5,0.5)/(1+sigmaq2))^2;
Lagrange3 = (integral(IntSx3,-0.5,0.5)/(1+sigmaq3))^2;

% Matched receive filters
Hfsq1 = @(f) sqrt((Lagrange1.*Sx(f))/sigmaq1)-Lagrange1;
Hfsq2 = @(f) sqrt((Lagrange2.*Sx(f))/sigmaq2)-Lagrange2;
Hfsq3 = @(f) sqrt((Lagrange3.*Sx(f))/sigmaq3)-Lagrange3;

figure(1);
plot(fuh,Hfsq1(fuh),fuh,Hfsq2(fuh),fuh,Hfsq3(fuh));
title('Frequency response of matched receive filters');
xlabel('Normalized frequency (f)');
ylabel('H(f)');
legend('0.75 bits/sample','2 bits/sample','5 bits/sample');

% Matched transmit filters
Gfsq1 = @(f) sqrt(sigmaq1./(Lagrange1.*Sx(f)))-(sigmaq1./Sx(f));
Gfsq2 = @(f) sqrt(sigmaq2./(Lagrange2.*Sx(f)))-(sigmaq2./Sx(f));
Gfsq3 = @(f) sqrt(sigmaq3./(Lagrange3.*Sx(f)))-(sigmaq3./Sx(f));

figure(2);
plot(fuh,Gfsq1(fuh),fuh,Gfsq2(fuh),fuh,Gfsq3(fuh));
title('Frequency response of matched transmit filters');
xlabel('Normalized frequency (f)');
ylabel('G(f)');
legend('0.75 bits/sample','2 bits/sample','5 bits/sample');

% Output noise PSD
Nf1 = sigmaq1.*Hfsq1(fuh);
Nf2 = sigmaq2.*Hfsq2(fuh);
Nf3 = sigmaq3.*Hfsq3(fuh);

% Output signal PSD
Yf1 = Sx(fuh).*Hfsq1(fuh).*Gfsq1(fuh);
Yf2 = Sx(fuh).*Hfsq2(fuh).*Gfsq2(fuh);
Yf3 = Sx(fuh).*Hfsq3(fuh).*Gfsq3(fuh);

figure(3);
semilogy(fuh,Yf1,'c');
hold on;
semilogy(fuh,Nf1);
semilogy(fuh,Yf2,'m');
semilogy(fuh,Nf2,'g');
semilogy(fuh,Yf3,'k');
semilogy(fuh,Nf3,'r');
title('Output signal and noise power spectrums');
xlabel('Normalized frequency (f)');
ylabel('Y(f), N(f) (logarithmic scale)');
legend('Signal (0.75 bits/sample)','Noise (0.75 bits/sample)','Signal (2 bits/sample)','Noise (2 bits/sample)','Signal (5 bits/sample)','Noise (5 bits/sample)');
hold off;

% Signal-to-Noise ratios
SNR1 = 10*log10(sum(Yf1)/sum(Nf1));
SNR2 = 10*log10(sum(Yf2)/sum(Nf2));
SNR3 = 10*log10(sum(Yf3)/sum(Nf3));

%------------ Exercise 2b ------------%

% Frequency sampling of the matched receive filter at low resolution
Fh1 = sqrt(Hfsq1(flf));
Fh2 = sqrt(Hfsq2(flf));
Fh3 = sqrt(Hfsq3(flf));
hn1 = FrSamp(Fh1);
hn2 = FrSamp(Fh2);
hn3 = FrSamp(Fh3);

% Frequency sampling of the matched transmit filter at low resolution
Fg1 = sqrt(Gfsq1(flf));
Fg2 = sqrt(Gfsq2(flf));
Fg3 = sqrt(Gfsq3(flf));
gn1 = FrSamp(Fg1);
gn2 = FrSamp(Fg2);
gn3 = FrSamp(Fg3);

figure(4);
hold on;
stem(0:1:NL-1,hn3,'r');
stem(0:1:NL-1,hn2,'g');
stem(0:1:NL-1,hn1,'b');
title(['Impulse response of matched receive filters (sampled from ' int2str(NL) ' points)']);
xlabel('n');
ylabel('h(n)');
legend('5 bits/sample','2 bits/sample','0.75 bits/sample');
axis([0 NL-1 0 0.6]);
hold off;

figure(5);
hold on;
stem(0:1:NL-1,gn3,'r');
stem(0:1:NL-1,gn2,'g');
stem(0:1:NL-1,gn1,'b');
title(['Impulse response of matched transmit filters (sampled from ' int2str(NL) ' points)']);
xlabel('n');
ylabel('g(n)');
legend('5 bits/sample','2 bits/sample','0.75 bits/sample');
axis([0 NL-1 -0.5 2.1]);
hold off;

% Frequency response of calculated receive filters
Hf1 = fft(hn1);
Hf2 = fft(hn2);
Hf3 = fft(hn3);
% High resolution frequency response of exact receive filters
Hfe1 = sqrt(Hfsq1(fhf(fhi)));
Hfe2 = sqrt(Hfsq2(fhf(fhi)));
Hfe3 = sqrt(Hfsq3(fhf(fhi)));

figure(6);
plot(flf(fli),abs(Hf1(fli)),fhf(fhi),abs(Hfe1),flf(fli),abs(Hf2(fli)), fhf(fhi),abs(Hfe2),flf(fli),abs(Hf3(fli)),fhf(fhi),abs(Hfe3));
title('Frequency response of exact and calculated receive filters');
xlabel('Normalized frequency (f)');
ylabel('H(f)');
legend('Calculated (0.75 bits/sample)','Exact (0.75 bits/sample)','Calculated (2 bits/sample)','Exact (2 bits/sample)','Calculated (5 bits/sample)','Exact (5 bits/sample)');

% Frequency response of calculated transmit filters
Gf1 = fft(gn1);
Gf2 = fft(gn2);
Gf3 = fft(gn3);
% High resolution frequency response of exact receive filters
Gfe1 = sqrt(Gfsq1(fhf(fhi)));
Gfe2 = sqrt(Gfsq2(fhf(fhi)));
Gfe3 = sqrt(Gfsq3(fhf(fhi)));

figure(7);
plot(flf(fli),abs(Gf1(fli)),fhf(fhi),abs(Gfe1),flf(fli),abs(Gf2(fli)), fhf(fhi),abs(Gfe2),flf(fli),abs(Gf3(fli)),fhf(fhi),abs(Gfe3));
title('Frequency response of exact and calculated transmit filters');
xlabel('Normalized frequency (f)');
ylabel('H(f)');
legend('Calculated (0.75 bits/sample)','Exact (0.75 bits/sample)', 'Calculated (2 bits/sample)','Exact (2 bits/sample)','Calculated (5 bits/sample)','Exact (5 bits/sample)');

% Signal modification
Sigmod1 = @(f) 1-sqrt(Lagrange1.*(sigmaq1./Sx(f)));
Sigmod2 = @(f) 1-sqrt(Lagrange2.*(sigmaq2./Sx(f)));
Sigmod3 = @(f) 1-sqrt(Lagrange3.*(sigmaq3./Sx(f)));

figure(8);
plot(fhf(fhi),abs(Sigmod1(fhf(fhi))),fhf(fhi),abs(Sigmod2(fhf(fhi))), fhf(fhi),abs(Sigmod3(fhf(fhi))));
title('Signal modification through the system');
xlabel('Normalized frequency (f)');
ylabel('|H(f)G(f)|');
legend('0.75 bits/sample','2 bits/sample','5 bits/sample');

%------------ Exercise 2c ------------%

% Generating the input signal
n = 100000;
xn = oppgave1c(n,rho);
un1 = conv(xn,gn1);
un2 = conv(xn,gn2);
un3 = conv(xn,gn3);

% Quantize the signal
quantlimit = 6;
quantstep1 = sqrt(12*sigmaq1);
extreme1 = floor(quantlimit/quantstep1)*quantstep1;
quantlevels1 = [-(extreme1+quantstep1):quantstep1/2:(extreme1+quantstep1)];
[index,quants1] = quantiz(un1,quantlevels1(2:2:length(quantlevels1)),quantlevels1(1:2:length(quantlevels1)));

quantstep2 = sqrt(12*sigmaq2);
extreme2 = floor(quantlimit/quantstep2)*quantstep2;
quantlevels2 = [-(extreme2+quantstep2):quantstep2/2:(extreme2+quantstep2)];
[index,quants2] = quantiz(un2,quantlevels2(2:2:length(quantlevels2)),quantlevels2(1:2:length(quantlevels2)));

quantstep3 = sqrt(12*sigmaq3);
extreme3 = floor(quantlimit/quantstep3)*quantstep3;
quantlevels3 = [-(extreme3+quantstep3):quantstep3/2:(extreme3+quantstep3)];
[index,quants3] = quantiz(un3,quantlevels3(2:2:length(quantlevels3)),quantlevels3(1:2:length(quantlevels3)));

figure(9);
subplot(1,3,1);
hist(quants1,quantlevels1(1:2:length(quantlevels1)));
subplot(1,3,2);
hist(quants2,quantlevels2(1:2:length(quantlevels2)));
subplot(1,3,3);
hist(quants3,quantlevels3(1:2:length(quantlevels3)));
annotation('textbox', [0 0.9 1 0.1],'String', 'Quantized output from transmit filters with 0.75, 2 and 5 bits/sample','EdgeColor', 'none','HorizontalAlignment', 'center');

% Calculate entropy
quant1 = hist(quants1,quantlevels1(1:2:length(quantlevels1)));
quant2 = hist(quants2,quantlevels2(1:2:length(quantlevels2)));
quant3 = hist(quants3,quantlevels3(1:2:length(quantlevels3)));

entr1 = 0;
entr2 = 0;
entr3 = 0;
for i = 1:length(quant1)
    Pi = quant1(i)/sum(quant1);
    if (Pi>0)
        entr1 = entr1 -(Pi*log2(Pi));
    end
end
for i = 1:length(quant2)
    Pi = quant2(i)/sum(quant2);
    if (Pi>0)
        entr2 = entr2 -(Pi*log2(Pi));
    end
end
for i = 1:length(quant3)
    Pi = quant3(i)/sum(quant3);
    if (Pi>0)
        entr3 = entr3 -(Pi*log2(Pi));
    end
end

% Output signals
yn1 = conv(quants1,hn1);
yn2 = conv(quants2,hn2);
yn3 = conv(quants3,hn3);

% Simulated Signal-Noise Ratios
delay = floor(NL/2);
SNR_sim1 = 10*log10(sum(xn.^2)/sum((xn-yn1((2*delay)+1:1:length(yn1)-(2*delay))).^2));
SNR_sim2 = 10*log10(sum(xn.^2)/sum((xn-yn2((2*delay)+1:1:length(yn2)-(2*delay))).^2));
SNR_sim3 = 10*log10(sum(xn.^2)/sum((xn-yn3((2*delay)+1:1:length(yn3)-(2*delay))).^2));