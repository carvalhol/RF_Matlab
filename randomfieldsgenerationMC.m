function y = randomfieldsgenerationMC( wN, N, corrMod, Nmc, t )
%
% wN : cut-off frequency

% size
Nt = size(t);
if Nt(1)>Nt(2); t = t'; end
Nt = max(Nt);

% frequency discretization
wk = (0:N-1)'/N*wN;

% random phases
phik = reshape( permute( repmat( 2*pi*rand(N,Nmc), [1 1 Nt] ), ...
                                                    [1 3 2] ), N, Nmc*Nt );


% spectrum 
Sk = correlationModel( wk, corrMod);

% generation of random field
y = reshape( 2* sqrt( Sk * wN / N ) * cos( repmat(wk*t,[1 Nmc]) + phik ), Nt, Nmc );

%==========================================================================
% SPECTRUM DEPENDING ON CORRELATION MODEL
function Sk = correlationModel( wk, model )
switch lower(model)
    case 'cardinalsine2'
        Sk = tripuls(1/pi*wk');
        
    case 'gaussian'
        Sk = sqrt(pi/10) * exp( - 1/(4*10)*wk'.^2) ;
        
    case 'exponential'
        Sk = 2 * 0.1 ./ ( (0.1)^2 +  wk'.^2 ) ;
        
        
    case 'powerlaw'
        error('not implemented yet')
        
    case 'vonkarman'
        error('not implemented yet')
        
    case 'fractal'
        error('not implemented yet')
        
    otherwise
        error('this correlation model does not exist')
        
end
