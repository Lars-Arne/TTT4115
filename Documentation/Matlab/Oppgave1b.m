clear all

f=[-1/2:0.01:1/2];
rho=0.9;
Sx=(1-rho^2)./(1+rho^2-2*rho*cos(2*pi*f));
plot(f,Sx)
title('Power spectral density')
ylabel('Sx')
xlabel('f')