%% PAM-4 Transmitter
T = 1;
a = 0.3;
tstep = 0.1;
t = -10:tstep:10-tstep;
g = sinc(t/T) .* (cos(a*pi*t/T) ./ (1 - (2*a*t/T).^2));
figure(1)
plot(t,g)

%%
ak = zeros(1,length(t));
am = randsrc(1,11,[-3,-1,1,3]);
k = (length(t)/2:(1/tstep):length(t));
ak(k) = am;
figure(2)
stem(t,ak)

%%
s = conv(ak, g);
figure(3)
hold on
t2 = -20:tstep:20-2*tstep;
plot(t2,s)
stem(t,ak)
hold off

%% rx
N = 0.5*max(s);
%n = N*gaussmf(t2, [0.5 5]);
n = N*rand(1, length(t2));
figure(4)
plot(t2, n)
%% 
r = s+n;
gn = sinc((-t)/T) .* (cos(a*pi*(-t)/T) ./ (1 - (2*a*(-t)/T).^2));
figure()
plot(t2,r)
h = conj(gn);
figure()
plot(t,h)
%% 

t3 = -30:tstep:30-3*tstep;
y = conv(r, gn);
figure()
plot(t3,y)

%% amostragem
for R = 1:length(t3)
    if (mod(R,100) == 0)
        
    end
end
