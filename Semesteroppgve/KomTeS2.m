clear all;

%---------Exercise 1a---------%

% Filter coefficients
k = tan([0.39, 0.41, 0.45, 0.49]*pi);

% blank input vectors
for i=[1:8]
    sz(i,:) = zeros(1,99);
end

% output signals with unit impulse input to one channel at a time
for i=[1:8]
    s = sz;
    s(i,:) = [1 zeros(1,98)];
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
linspecs = ['--';': ';'--';': ';'--';': ';'--';': '];

hFig = figure(1);
set(hFig, 'Position', [100 100 1000 600]);
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
set(hFig, 'Position', [100 100 1000 600]);
plot(frscale,abs(C(1,:)),linspecs(1,:),frscale,abs(C(2,:)),linspecs(2,:),frscale,abs(C(3,:)),linspecs(3,:),frscale,abs(C(4,:)),linspecs(4,:),frscale,abs(C(5,:)),linspecs(5,:),frscale,abs(C(6,:)),linspecs(6,:),frscale,abs(C(7,:)),linspecs(7,:),frscale,abs(C(8,:)),linspecs(8,:));
legend('CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8');
xlabel('Normalized frequency (f)');
ylabel('Frequency response X(f)');
title('Frequency response per channel after transmission, before filtering');

%---------Exercise 2b---------%

[u1,u2,u3,u4,u5,u6,u7,u8] = RFB(c(8,:),k);