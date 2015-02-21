clear all


n=10000;
tau=[-50:50];
rho=0.9;
b=sqrt(1-rho^2);
a=[1 -rho];
e=randn(1,n);
x=filter(b,a,e);
acorr=xcorr(x);



Sx=fftshift(fft(acorr)./n); %What is this?



figure(1)
plot(tau,acorr(n-50:n+50)./n,tau,rho.^(abs(tau)));
title('Oppgave1d')
ylabel('Rx(tau)')
xlabel('tau')
figure(2)
plot([-n+1:n-1],abs(Sx))
