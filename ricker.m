% Ricker Wavelet generation, in time & frequency domain 
% Parameters: 
% f0 = central (carrier-fundamental) frequency
% t0 = delay (the position of the pick in time domain)

function f = ricker(t,f0,t0)
arg = pi*f0*(t-t0);
arg = arg.*arg;
f = (1-2*arg).*exp(-arg);


% % Ricker wavelet in frequency domain
% function f = ricker(w,f0,t0)
% 
% f = ( 2*sqrt(pi)*exp(-1i*t0.*w)/(pi*f0) ) .* ( (w/(2*pi*f0)).*(w/(2*pi*f0)) )...
%     .* exp(-(w/(2*pi*f0)).*(w/(2*pi*f0)) );
% end
