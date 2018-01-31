function out = mkrf1D(method, PSD, x, lc, seed)
    %method :
    % 1-spectral with contant amplitude
    % 2-spectral with random amplitude
    % 3-randomization
    %PSD - Power Spectral Density: 'gaussian'
    %x - vector of spatial coordinates
    %lc - double of correlation length
    %seed - random seed

    switch method

        case 1
            out = rf1d_structured(PSD, x, lc, seed);
        case 2
            out = rf1d_structured(PSD, x, lc, seed, true);
        case 3
            out = rf1d_randomization(PSD, x, lc, seed);
        otherwise
            error('Options are: 1-spectral with contant amplitude, 2-spectral with random amplitude 3-randomization.')
    end    
    

end


function out = rf1d_structured(PSD, x, lc, seed, ampRand)
    % constants
    xrange = max(x)-min(x);
    Nfft = 8;
    
    if nargin<5; ampRand = false; end
    
    % discretization in fourier space
    dk = min( [pi/2/lc ; pi/xrange] )*.95;
    kmax = 4*pi/lc+dk;
    k = (0:dk:kmax);
    Nk = length(k);
    xk = linspace( 0, 2*pi/dk, Nfft*Nk );

    % construction of power spectral density
    switch PSD

        case 'gaussian'
            Sk = (lc * exp(-((k.^2) * (lc^2))/(4*pi)))';
        
        % triangular power spectrum = sinc^2 correlation
        case 'sinc2'
            Sk = tripuls( k, 4*pi/Lc )';
            if ~isempty(dL) && (dL>Lc)
                Sk ( k>2*pi/dL, : ) = 0;
            end
        otherwise
            error('Options are: gaussian. Other first-order marginal law have not been implemented yet')
    end

    Nmc = 1;

    % seeds in space
    rng( seed );
    phi = randn( Nk, Nmc );
    %out = 2*pi*rand( Nk, Nmc );

    % random field in fourier space (isotropic case)
    if(ampRand)
        ksi = randn( Nk, Nmc );
        Sk_mat = ksi .* repmat(sqrt(Sk), [1 Nmc]);
    else
        Sk_mat = repmat(sqrt(Sk), [1 Nmc]);
    end
    cosPhi = fft(phi, [], 1);
    cosPhi = cosPhi .* Sk_mat;

    % return to real space
    out = ifft(cosPhi, Nfft*Nk, 1, 'symmetric' );
    indx = find(xk>xrange, 1);
    

    amp = Nfft*Nk*dk*lc/2/pi/sqrt(2);
    out = amp * interp1(xk(1:indx), out(1:indx,:), x-min(x));
    %======= MAKEGAUSS ========
    % generator of gaussian random field with unit variance and zero mean
    % R. Cottereau 08/2012
end

function out = rf1d_randomization(PSD, x, lc, seed)
    % construction of the vector k
    switch PSD

        case 'gaussian'
            ksi = randn( Nk, Nmc );
        
        otherwise
            error('Options are: gaussian. Other first-order marginal law have not been implemented yet')
    end
    
    dk = min( [pi/2/lc ; pi/xrange] )*.95;
    kmax = 4*pi/lc+dk;
    k = (0:dk:kmax);
    Nk = length(k);
    xk = linspace( 0, 2*pi/dk, Nfft*Nk );

end