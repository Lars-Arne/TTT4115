function [u1,u2,u3,u4,u5,u6,u7,u8]=RFB(y,k)
%Receiver filter bank which splits and demodulates the signal vector y into 
%the vectors u1,u2,...,u8  
%k is the vector of 4 filter coefficients

[Yu5,Yu6]=anaN(y,k );

[Yu1,Yu2]=anaN(Yu5,k);
[Yu3,Yu4]=anaN(Yu6,k);

[xu1,xu2]=anaN(Yu1,k );
[xu3,xu4]=anaN(Yu2,k );
[xu5,xu6]=anaN(Yu3,k );
[xu7,xu8]=anaN(Yu4,k );
 
%%Aligning with inputs (Delay=6)
u1=xu1(8:length(xu1)-7);
u2=xu2(8:length(xu2)-7);
u3=xu3(8:length(xu1)-7);
u4=xu4(8:length(xu2)-7);

u5=xu5(8:length(xu1)-7);
u6=xu6(8:length(xu2)-7);
u7=xu7(8:length(xu1)-7);
u8=xu8(8:length(xu2)-7);