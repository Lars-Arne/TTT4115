function y=synN(x1,x2,k)
 
%Synthesis filter block with four lattice stages
%Input signals: x1,x2
%Filter coefficients contained in array k
%Output signal: y
 
[s1,s2]=syn(x1,x2,k(4));
[v1,v2]=syn(s1,s2,k(3));
[v3,v4]=syn(v1,v2,k(2));
[v5,v6]=syn(v3,v4,k(1));
y=(upsample(v5,2)+upsample(v6,2,1));