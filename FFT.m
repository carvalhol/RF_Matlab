A = 1
x = 1:0.01:10;
x = cos(A*x);

Fx = fft(x)

plot(Fx)