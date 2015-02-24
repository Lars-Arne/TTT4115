clear all

n=0:0.01:10000;

A=randn(1,length(n));

fc=1;
X=sin(2*pi*fc*n);

Y=sum(A.*X);