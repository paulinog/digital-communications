close all;
clc; clearvars;
%%
f0 = 25;
fc = 1e3;
fs = 8e3;

t = 0:1/fs:0.1;

x = sin(2*pi*f0*t);
xrf = x.*cos(2*pi*fc*t);

[num, den] = butter(10,fc*2/fs);
y = amdemod(xrf, fc, fs,0,0,num, den);

figure()
hold on
plot(t,x)
plot(t,xrf)
plot(t, y, '--')
legend('x','x_{RF}','y')
title('AM Demod')

%%


