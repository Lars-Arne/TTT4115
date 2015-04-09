clear all;
graphPos = [100 100 1000 600];
nz = 99;
sone = [1 zeros(1,nz-1)];
nscale = [0:0.1:0.1*(nz-1)];

% blank input vectors
for i=[1:8]
    sz(i,:) = zeros(1,nz);
end

linspecs = ['--';': ';'--';': ';'--';': ';'--';': '];

%---------Exercise 1a---------%

% Filter coefficients
k = tan([0.39, 0.41, 0.45, 0.49]*pi);

% output signals with unit impulse input to one channel at a time
for i=[1:8]
    s = sz;
    s(i,:) = sone;
    x(i,:) = TFB(s(1,:),s(2,:),s(3,:),s(4,:),s(5,:),s(6,:),s(7,:),s(8,:),k);
end

% Check orthogonality through dot products between all combinations of
% inputs. signals are orthogonal if their dot product is 0.
% 'dots' matrix should have 1s along diagonal and 0 everywhere else if
% signals are orthogonal.
for i=[1:8]
    for j=[1:8]
        dots(i,j) = dot(x(i,:),x(j,:));
    end
end

dots

%---------Exercise 1b---------%

frscale = 0:1/length(x(1,:)):1-1/length(x(1,:));

hFig = figure(1);
set(hFig, 'Position', graphPos);
for i=[1:8]
    freq(i,:) = fftshift(fft(x(i,:)));
end
plot(frscale,abs(freq(1,:)),linspecs(1,:),frscale,abs(freq(2,:)),linspecs(2,:),frscale,abs(freq(3,:)),linspecs(3,:),frscale,abs(freq(4,:)),linspecs(4,:),frscale,abs(freq(5,:)),linspecs(5,:),frscale,abs(freq(6,:)),linspecs(6,:),frscale,abs(freq(7,:)),linspecs(7,:),frscale,abs(freq(8,:)),linspecs(8,:));
legend('CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8');
xlabel('Normalized frequency (f)');
ylabel('Frequency response X(f)');
title('Frequency response per channel before transmission');

%---------Exercise 2a---------%

h = [0.3336,0.2975,0.1328,0.0729,0.0389,0.0274,0.0172,0.0140,0.0098,0.0087,0.0064,0.0061,0.0047,0.0048,0.0037,0.0042,0.0029,0.0046,0.0010,0.0086];

for i=[1:8]
    c(i,:) = filter(h,1,x(i,:));
    C(i,:) = fftshift(fft(c(i,:)));
end

hFig = figure(2);
set(hFig, 'Position', graphPos);
plot(frscale,abs(C(1,:)),linspecs(1,:),frscale,abs(C(2,:)),linspecs(2,:),frscale,abs(C(3,:)),linspecs(3,:),frscale,abs(C(4,:)),linspecs(4,:),frscale,abs(C(5,:)),linspecs(5,:),frscale,abs(C(6,:)),linspecs(6,:),frscale,abs(C(7,:)),linspecs(7,:),frscale,abs(C(8,:)),linspecs(8,:));
legend('CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8');
xlabel('Normalized frequency (f)');
ylabel('Frequency response X(f)');
title('Frequency response per channel after transmission, before equalization');

%---------Exercise 2b---------%

[u1a,u2a,u3a,u4a,u5a,u6a,u7a,u8a] = RFB(c(4,:),k);
[u1b,u2b,u3b,u4b,u5b,u6b,u7b,u8b] = RFB(c(8,:),k);

hFig = figure(3);
set(hFig, 'Position', graphPos);
u8n = stem(u8b(1:20));
hold on;
u4n = stem(u4a(1:20),'r');
hold off;
legend([u8n,u4n],'CH8','CH4');
xlabel('n');
ylabel('u[n]');
title('Impulse responses for channels 4 and 8 after transmission, before equalization');

%---------Exercise 2c---------%

testvect = besselj(0,nscale);
xtest = TFB(sz(1,:),sz(2,:),sz(3,:),testvect,sz(5,:),sz(6,:),sz(7,:),sz(8,:),k);
ctest = filter(h,1,xtest);
vtest = filter(1,h,ctest);
[u1,u2,u3,u4,u5,u6,u7,u8] = RFB(vtest,k);

noise = sum(testvect-u4)

hFig = figure(4);
set(hFig, 'Position', graphPos);
zplane(h,1);
title('Pole-zero plot of channel');

%---------Exercise 2d---------%

for i=[1:8]
    v(i,:) = filter(1,h,c(i,:));
    V(i,:) = fftshift(fft(v(i,:)));
end

[u1a,u2a,u3a,u4a,u5a,u6a,u7a,u8a] = RFB(v(4,:),k);
[u1b,u2b,u3b,u4b,u5b,u6b,u7b,u8b] = RFB(v(8,:),k);

hFig = figure(5);
set(hFig, 'Position', graphPos);
u8n = stem(u8b(1:20));
hold on;
u4n = stem(u4a(1:20),':r');
hold off;
legend([u8n,u4n],'CH8','CH4');
xlabel('n');
ylabel('u[n]');
title('Impulse responses for channels 4 and 8 after transmission and equalization');

hFig = figure(6);
set(hFig, 'Position', graphPos);
plot(frscale,abs(V(1,:)),linspecs(1,:),frscale,abs(V(2,:)),linspecs(2,:),frscale,abs(V(3,:)),linspecs(3,:),frscale,abs(V(4,:)),linspecs(4,:),frscale,abs(V(5,:)),linspecs(5,:),frscale,abs(V(6,:)),linspecs(6,:),frscale,abs(V(7,:)),linspecs(7,:),frscale,abs(V(8,:)),linspecs(8,:));
legend('CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8');
xlabel('Normalized frequency (f)');
ylabel('Frequency response V(f)');
title('Frequency response per channel after transmission and equalization');

%---------Exercise 2e---------%

for i=[1:8]
    w = 0.02*randn(1,length(c(i,:)));
    v(i,:) = filter(1,h,c(i,:)+w);
    %V(i,:) = fftshift(fft(v(i,:)));
    [u(1,:),u(2,:),u(3,:),u(4,:),u(5,:),u(6,:),u(7,:),u(8,:)] = RFB(v(i,:),k);
    sr(i,:) = u(i,:);
    sigmaq(i) = sum(abs(sr(i,:)-sone));
end

sigmaq

%---------Exercise 2f---------%

% Find center frequency of channel 3
[Y,I] = max(abs(freq(3,1:length(freq(3,:))/2)));
cf3 = frscale(I);

% Generate cosine with center frequency of channel 3
cosvect = cos(2*pi*cf3.*nscale);

% Transmit cosine in channel 3, zero in other channels
s = sz;
s(3,:) = cosvect;
xcos = TFB(s(1,:),s(2,:),s(3,:),s(4,:),s(5,:),s(6,:),s(7,:),s(8,:),k);
ccos = filter(h,1,xcos);
[u1c,u2c,u3c,u4c,u5c,u6c,u7c,u8c] = RFB(ccos,k);