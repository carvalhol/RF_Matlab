function y = randomgenerationfieldsFFTMeth( t,corrMod, L )
%
% wN : cut-off frequency
Nt = size(t);
if Nt(1)>Nt(2); t = t'; end
Nt = max(Nt);
T = max(t)-min(t);
dw = 2*pi/T;
w = (0:Nt-1)*dw;

% random phases
phik = rand(1,Nt) ; 

% spectrum
Sk = correlationModel( w, corrMod, L );
Dk = sqrt( dw*Sk ).*exp(complex(0,1)*2*pi*phik);
Dk = [Dk conj(Dk(end:-1:1))];

%Computes the inverse of fft
y = real((2*Nt)*ifft(Dk)) ; 
y = y(1:Nt);

%==========================================================================
% SPECTRUM DEPENDING ON CORRELATION MODEL
function Sk = correlationModel( wk, model,L )
switch lower(model)
    case 'cardinalsine2'
        eta = pi*L ; 
        Sk = 1/(eta*pi)*tripuls( wk,2/eta) ; 
        
    case 'gaussian'
        eta = L/sqrt(pi); 
        Sk = 1/(2*pi)*exp(-eta^2/4*wk.^2);
        
    case 'exponential'
        Sk = 1/pi*(L/2)./(1+(L/2*wk).^2);
        
    otherwise
        error('this correlation model does not exist')
        
end
