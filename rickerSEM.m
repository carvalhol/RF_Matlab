% Ricker Wavelet generation, in time & frequency domain 
% Parameters: 
% tau = central (carrier-fundamental) frequency
% t0 = delay (the position of the pick in time domain)

function signal = rickerSEM(time,f0,tau, amplitude)
alpha = -1*pi^2*f0^2;
sigma = alpha * (time - tau).^2;
signal = amplitude*(time < 2.5*tau)...
         .*2*alpha.*(1+2*sigma).*exp(sigma);


% % Ricker wavelet in frequency domain
% function f = ricker(w,f0,t0)
% 
% f = ( 2*sqrt(pi)*exp(-1i*t0.*w)/(pi*f0) ) .* ( (w/(2*pi*f0)).*(w/(2*pi*f0)) )...
%     .* exp(-(w/(2*pi*f0)).*(w/(2*pi*f0)) );
% end
