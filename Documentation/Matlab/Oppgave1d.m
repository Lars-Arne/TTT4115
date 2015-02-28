clear all

n=100;
rho=0.9;
tau=[-50:50];


x=Oppgave1c(n,rho);
Rxsim=xcorr(x)./n;
Rxcalc=rho.^(abs(tau));
plot(tau,Rxsim(n-50:n+50),tau,Rxcalc)
title('Comparison of calculated and simulated autocorrelation of X')
ylabel('Rx')
xlabel('tau')
legend('Simulated','Calculated')