clear all
deltau=1;
H=0.75;
%nabla2=deltau^2/(log((2^(2*H))/(2*pi)));

nabla2=(2*pi*exp(1)*deltau)/2^(2*H);
sigmaq2=nabla2/12;



Gf=conj(Hf)*Sx./(abs(Hf).^2*Sx+sigmaq2);
