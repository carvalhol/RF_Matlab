clear all
clc

pad = 2^10;

t = linspace(0,.1,pad);

fa = 1/mean(diff(t));

w = 0:fa/pad:fa-fa/pad;

vp=380;
vs=150;

f0 = 380;
t0 = 1/100;

arg = pi*f0*(t-t0);
arg = arg.*arg;
f = (1-2*arg).*exp(-arg);
F = fft(f);
%%
figure(1);clf
plot(t,f)
xlabel('Time [s]')
ylabel('Amplitude [a.u]')
title('Ricker Pulse')
figure(2);clf
plot(w(1:end/2),abs(F(1:end/2)/max(F)))
xlim([0 w(end/2)])
xlabel('Frequency [Hz]')
ylabel('Amplitude [a.u]')
title('Ricker Pulse')

[~, b]=max(abs(fft(f)));

w(b)
lambda_p = vp/w(b)
lambda_s = vs/w(b)
