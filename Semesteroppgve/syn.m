function [v1,v2]=syn(u1,u2,k)
%Synthesis filter block
%Input signals: u1,u2
%Filter coefficient: k
%Output signals: v1,v2
 
v1=([0,u1]-k*[0,u2])/sqrt(1+k^2);
v2=(k*[u1,0]+[u2,0])/sqrt(1+k^2);
 