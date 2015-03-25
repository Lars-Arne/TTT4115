function [x1,x2]=anaN(y,k)
%Analysis filter block with four lattice stages
%Input signal: y
%Filter coeffcients contained in array k
%Output signals: x1,x2

u1=downsample(y,2,0); 
u2=downsample(y,2,1); 
[u3,u4]=ana(u1,u2,k(1));
[u5,u6]=ana(u3,u4,k(2));
[u7,u8]=ana(u5,u6,k(3));
[x1,x2]=ana(u7,u8,k(4));
 