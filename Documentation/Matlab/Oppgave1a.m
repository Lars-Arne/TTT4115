clear all

tau=[-100:100];
rho=0.9;

Rx=rho.^(abs(tau));
plot(tau,Rx)
title('Autocorrelation of x')
xlabel('Tau')
ylabel('Rx(Tau)')