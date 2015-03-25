function y=TFB(x1,x2,x3,x4,x5,x6,x7,x8,k)
%Transmitter filter bank which combines the signal vectors x1,x2,...,x8 to 
%the signal y by frequency multiplex
%k is the vector of 4 filter coefficients

Y1=synN(x1,x2,k);
Y2=synN(x3,x4,k);
Y3=synN(x5,x6,k);
Y4=synN(x7,x8,k);

Y5=synN(Y1,Y2,k);
Y6=synN(Y3,Y4,k);

y=synN(Y5,Y6,k);