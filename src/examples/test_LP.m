close all;
clc; clearvars;
%%
f0 = 1000;
f1 = 250;
fs = 80e3;
fc = 10e3;

t = 0:1/fs:0.1;
% x = [1 2]*cos(2*pi*[f0 f1]'.*t) + randn(size(t))/10;
x = sin(2*pi*f0*t);%.*sin(2*pi*fc*t);% + randn(size(t))/10;
s = ammod(x,fc,fs);
% y = lowpass(x,150,fs);
% [num,den] = butter(10,f1*2*pi/2);
[num,den] = butter(10,fc*2/fs);
y = amdemod(x,fc,fs,0,0,num,den);

figure()
hold on
plot(t,x)
plot(t,y)

%%
% fc = 256;
% fs = 8192;
% t = 0:1/fs:1;
% x = sin(2*pi*fc*t) + randn(size(t))/10;
% % t_norm = t/length(t);
% x = sin(2*pi*fc*t);
% y = lowpass(x,fc/32,fs);
% % y = lowpass(x,0.001);
% figure()
% hold on
% plot(t,x)
% plot(t,y)

