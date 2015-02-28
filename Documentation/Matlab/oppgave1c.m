function [ x ] = Oppgave1c( n, rho )
%[ x ] = Oppgave1c( n, rho )
%Generates an AR-function from the filter H(z)=(1-rho^2)/(1-rho*z^-1) of
%length n.

b=sqrt(1-rho^2);
a=[1 -rho];
e=randn(1,n);
x=filter(b,a,e);

end

