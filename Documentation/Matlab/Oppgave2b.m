clear all;

rho = 0.9;

% Low resolution frequency vector
NL = 29;
f = 0:1/NL:1-1/NL;
% High resolution frequency vector
NH = 999;
f2 = 0:1/NH:1-1/NH;

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
Sx = @(f) 0.19./(1.81-(1.8.*cos(2*pi.*f)));

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

% Matched transmit filters
Gfsq1 = @(f) sqrt(sigmaq1./(Lagrange1.*Sx(f)))-(sigmaq1./Sx(f));
Gfsq2 = @(f) sqrt(sigmaq2./(Lagrange2.*Sx(f)))-(sigmaq2./Sx(f));
Gfsq3 = @(f) sqrt(sigmaq3./(Lagrange3.*Sx(f)))-(sigmaq3./Sx(f));

% Frequency sampling of the matched receive filter at low resolution
Fh1 = sqrt(Hfsq1(f));
Fh2 = sqrt(Hfsq2(f));
Fh3 = sqrt(Hfsq3(f));
hn1 = FrSamp(Fh1);
hn2 = FrSamp(Fh2);
hn3 = FrSamp(Fh3);

% Frequency sampling of the matched transmit filter at low resolution
Fg1 = sqrt(Gfsq1(f));
Fg2 = sqrt(Gfsq2(f));
Fg3 = sqrt(Gfsq3(f));
gn1 = FrSamp(Fg1);
gn2 = FrSamp(Fg2);
gn3 = FrSamp(Fg3);

figure(1);
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

figure(2);
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

% Half spectrum x-axis vector
fl = 1:ceil(length(f)/2);
fh = 1:ceil(length(f2)/2);

% Frequency response of calculated receive filters
Hf1 = fft(hn1);
Hf2 = fft(hn2);
Hf3 = fft(hn3);
% High resolution frequency response of exact receive filters
Hfe1 = sqrt(Hfsq1(f2(fh)));
Hfe2 = sqrt(Hfsq2(f2(fh)));
Hfe3 = sqrt(Hfsq3(f2(fh)));

figure(3);
plot(f(fl),abs(Hf1(fl)),f2(fh),abs(Hfe1),f(fl),abs(Hf2(fl)),f2(fh) ,abs(Hfe2),f(fl),abs(Hf3(fl)),f2(fh),abs(Hfe3));
title('Frequency response of exact and calculated receive filters');
xlabel('Normalized frequency (f)');
ylabel('H(f)');
legend('Calculated (0.75 bits/sample)','Exact (0.75 bits/sample)','Calculated (2 bits/sample)','Exact (2 bits/sample)','Calculated (5 bits/sample)','Exact (5 bits/sample)');

% Frequency response of calculated transmit filters
Gf1 = fft(gn1);
Gf2 = fft(gn2);
Gf3 = fft(gn3);
% High resolution frequency response of exact receive filters
Gfe1 = sqrt(Gfsq1(f2(fh)));
Gfe2 = sqrt(Gfsq2(f2(fh)));
Gfe3 = sqrt(Gfsq3(f2(fh)));

figure(4);
plot(f(fl),abs(Gf1(fl)),f2(fh),abs(Gfe1),f(fl),abs(Gf2(fl)),f2(fh), abs(Gfe2),f(fl),abs(Gf3(fl)),f2(fh),abs(Gfe3));
title('Frequency response of exact and calculated transmit filters');
xlabel('Normalized frequency (f)');
ylabel('H(f)');
legend('Calculated (0.75 bits/sample)','Exact (0.75 bits/sample)','Calculated (2 bits/sample)','Exact (2 bits/sample)','Calculated (5 bits/sample)','Exact (5 bits/sample)');

% Signal modification
Sigmod1 = @(f) 1-sqrt(Lagrange1.*(sigmaq1./Sx(f)));
Sigmod2 = @(f) 1-sqrt(Lagrange2.*(sigmaq2./Sx(f)));
Sigmod3 = @(f) 1-sqrt(Lagrange3.*(sigmaq3./Sx(f)));
figure(5);
plot(f2(fh),abs(Sigmod1(f2(fh))),f2(fh),abs(Sigmod2(f2(fh))), f2(fh),abs(Sigmod3(f2(fh))));
title('Signal modification through the system');
xlabel('Normalized frequency (f)');
ylabel('|H(f)G(f)|');
legend('0.75 bits/sample','2 bits/sample','5 bits/sample');