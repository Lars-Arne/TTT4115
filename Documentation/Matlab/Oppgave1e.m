clear all
n=100;
rho=0.9;

f=[-1/2:1/100:1/2];
x=Oppgave1c(n,rho);
acorr=xcorr(x);
Sxsim=fftshift(fft(acorr(n-50:n+50))./n);
Sxcalc=(1-rho^2)./(1+rho^2-2*rho*cos(2*pi*f));

plot(f,abs(Sxsim),f,Sxcalc)
title('Comparison of calculated and simulated power spectral densities of X')
ylabel('Sx')
xlabel('f')
legend('Simulated','Calculated')
