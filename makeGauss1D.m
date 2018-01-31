function gx = makeGauss1D( type, x, Lc, dL, Nmc, sd )

% constants
xmax = max(x)-min(x);
Nfft = 8;

% discretization in fourier space
dk = min( [pi/2/Lc ; pi/xmax] )*.95;
kmax = 4*pi/Lc+dk;
k = (0:dk:kmax);
Nk = length(k);
xk = linspace( 0, 2*pi/dk, Nfft*Nk );

% construction of power spectral density
switch type
    
    % triangle power spectrum = sinc^2 correlation
    case 'sinc2'
        Sk = tripuls( k, 4*pi/Lc )';
        if ~isempty(dL) && (dL>Lc)
            Sk ( k>2*pi/dL, : ) = 0;
        end
        
    % granular power spectrum
    case 'granular'
        error('granular should be checked')
%         r = (0:(N-1))'/N * L/Lc;
%         C = CorrelationDensityTheoretical( .57, r );
%         Sk = fft( C )*L/N;        

    % not implemented yet
    otherwise
        error('this first-order marginal law has not been implemented yet')
end

% seeds in space
rng( sd );
gx = randn( Nk, Nmc );

% random field in fourier space (isotropic case)
gk = fft( gx, [], 1 );
gk = gk .* repmat( sqrt(Sk), [1 Nmc] );

% return to real space
amp = Nfft*Nk*dk*Lc/2/pi/sqrt(2);
gx = ifft( gk, Nfft*Nk, 1, 'symmetric' );
indx = find( xk>xmax, 1 );
gx = amp * interp1( xk(1:indx), gx(1:indx,:), x-min(x) );

%======= MAKEGAUSS ========
% generator of gaussian random field with unit variance and zero mean
% R. Cottereau 08/2012