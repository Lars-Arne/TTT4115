function [v1,v2]=ana(u1,u2,k)
%Analysis filter block
%Input signals: u1,u2
%Filter coefficient: k
%Output signals: v1,v2
 
v1=(k*[0,u2]+[u1,0])/sqrt(1+k^2);
v2=([0,u2]-k*[u1,0])/sqrt(1+k^2);
 
 