clear all;
close all;

% Graph and plot settings
graphPos = [100 100 1000 600];
mpcols = 2;
mprows = 4;
linspecs = ['--';': ';'--';': ';'--';': ';'--';': '];

% Time domain vectors
nz = 1000;
sone = [1 zeros(1,nz-1)];
nscale = [0:0.01:0.01*(nz-1)];

% blank input vectors
sz = zeros(8,nz);

% Simulated channel and equalization filters
h = [0.3336,0.2975,0.1328,0.0729,0.0389,0.0274,0.0172,0.0140,0.0098,0.0087,0.0064,0.0061,0.0047,0.0048,0.0037,0.0042,0.0029,0.0046,0.0010,0.0086];
csim = dfilt.df1(h,1);
cinv = dfilt.df1(1,h);

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

for i=[1:8]
    freq(i,:) = fft(x(i,:));
end
hFig = figure(1);
set(hFig, 'Position', graphPos);
plot(frscale,abs(freq(1,:)),linspecs(1,:),frscale,abs(freq(2,:)),linspecs(2,:),frscale,abs(freq(3,:)),linspecs(3,:),frscale,abs(freq(4,:)),linspecs(4,:),frscale,abs(freq(5,:)),linspecs(5,:),frscale,abs(freq(6,:)),linspecs(6,:),frscale,abs(freq(7,:)),linspecs(7,:),frscale,abs(freq(8,:)),linspecs(8,:));
legend('CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8');
xlabel('Normalized frequency (f)');
ylabel('Frequency response X(f)');
title('Frequency response per channel before transmission');

%---------Exercise 2a---------%

for i=[1:8]
    c(i,:) = filter(csim,x(i,:));
    C(i,:) = fft(c(i,:));
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
ctest = filter(csim,xtest);
vtest = filter(cinv,ctest);
[u1,u2,u3,u4,u5,u6,u7,u8] = RFB(vtest,k);

noise = sum(testvect-u4)

hFig = figure(4);
set(hFig, 'Position', graphPos);
zplane(h,1);
title('Pole-zero plot of channel');

%---------Exercise 2d---------%

for i=[1:8]
    v(i,:) = filter(cinv,c(i,:));
    V(i,:) = fft(v(i,:));
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

% Find center frequencies
for i=[1:8]
    [Y,I] = max(abs(freq(i,1:length(freq(i,:))/2)));
    cf(i) = frscale(I);
end

% Channel to send sine on
testch = 4;

% Transmit sinusoid in select channel, zero in other channels
s = sz;
s(testch,:) = sin(2*pi*cf(testch)*nscale);
xcos = TFB(s(1,:),s(2,:),s(3,:),s(4,:),s(5,:),s(6,:),s(7,:),s(8,:),k);
ccos = filter(h,1,xcos);
[uc(1,:),uc(2,:),uc(3,:),uc(4,:),uc(5,:),uc(6,:),uc(7,:),uc(8,:)] = RFB(ccos,k);

Ccos = fft(xcos);

hFig = figure(8);
set(hFig, 'Position', graphPos);
for i=[1:8]
    subplot(mprows,mpcols,i);
    plot(nscale,uc(i,:));
    title(['Channel ' num2str(i)]);
end
%plot(nscale,sinvect,linspecs(1,:),nscale,u4c,linspecs(2,:));

%---------Exercise 3a---------%

Cf = exp(-4.6.*abs(cf));
%mintheta = (1/(prod(abs(Cf).^2)))^(1/8);
mintheta = 1;
theta = mintheta:10*mintheta:(1000*mintheta)-mintheta;
%for i=1:length(theta)
%    Cap(i) = 0.5*sum(log2(theta(i).*(abs(Cf).^2)));
%end

for i=1:8
    Cap(i,:) = 0.5.*log2((abs(Cf(i))^2).*theta);
end

Levels = round(2.^Cap);

hFig = figure(9);
set(hFig, 'Position', graphPos);
%plot(theta,Levels(1,:),linspecs(1,:),theta,Levels(2,:),linspecs(2,:),theta,Levels(4,:),linspecs(3,:),theta,Levels(4,:),linspecs(4,:),theta,Levels(5,:),linspecs(5,:),theta,Levels(6,:),linspecs(6,:),theta,Levels(7,:),linspecs(7,:),theta,Levels(8,:),linspecs(8,:));
%legend('CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8','Location','northwest');
%grid;
bar3(flipud(Levels),1);
set(gca,'YTickLabel',{'8','7','6','5','4','3','2','1'});
set(gca,'PlotBoxAspectRatio',[2 1 1]);
set(gca,'XTick',[1:10:length(theta)]);
set(gca,'XTickLabel',theta(1:10:end));
title('Optimal number of levels per channel');
ylabel('Channel');
xlabel('theta');
zlabel('Levels');

%---------Exercise 3b---------%

for i=1:8
    Sx(i,:) = theta-(1/(abs(Cf(i))^2));
end

Sx(Sx<0) = 0;

hFig = figure(10);
set(hFig, 'Position', graphPos);
bar3(flipud(Sx),1);
set(gca,'YTickLabel',{'8','7','6','5','4','3','2','1'});
set(gca,'PlotBoxAspectRatio',[2 1 1]);
set(gca,'XTick',[1:10:length(theta)]);
set(gca,'XTickLabel',theta(1:10:end));
title('Optimal power spectral density per channel');
ylabel('Channel');
xlabel('theta');
zlabel('Sx(f)');

%---------Exercise 3c---------%

theta3index = find(Levels(4,:) == 3,1,'first');
theta3level = theta(theta3index);

% Generate gaussian input with unit variance in channels 4 and 8
s = sz;
s(4,:) = round(30.*randn(1,nz));
s(8,:) = round(5.*randn(1,nz));

x = TFB(s(1,:),s(2,:),s(3,:),s(4,:),s(5,:),s(6,:),s(7,:),s(8,:),k);
c = filter(csim,x);
C = fft(c);

hFig = figure(11);
set(hFig, 'Position', graphPos);
plot(frscale,abs(C));
xlabel('Normalized frequency (f)');
ylabel('Power spectral density Sx(F)');
title('Power spectral density of transmitted signal');

y = filter(cinv,c);
[u1,u2,u3,u4,u5,u6,u7,u8] = RFB(x,k);

% Calculate bit error rate
ber4 = sum(abs(s(4,:)-u4)>0.5)
ber8 = sum(abs(s(8,:)-u8)>0.5)