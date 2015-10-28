function k = randomField( law, correl, L, mu, s, Nmc, x, y, z, dL, sd, change)
% RANDOMFIELD to generate realizations of a random field with given
% first-order marginal law and correlation structure
%
% syntax: E = randomField( law, correlation, L, mu, s, N, x,[y],[z],dL,sd)
%
%  law: 1st marginal law of the process: 'gamma', 'uniform', 'gauss'
%       and 'lognormal' implemented only
%  correlation: type of correlation model: 'sinc2' (=triangle spectrum),
%               and 'granular' implemented only
%  L: correlation lengths [1*d real]. If L is a single scalar, the same
%     value is used in all directions
%  mu: mean [scalar]
%  s: variance [scalar]
%  N: number of realizations to be computed [integer]
%  x,y,z: nodal coordinates in each direction [vectors]
%  dL: smoothing length (when dL>L, the smoothing length is considered to
%      upscale the random field at dL)
%  sd: integer to initiate the sequence of quasi-random numbers used to
%      generate the random fields (see function RNG). By default, the
%      generate is set to rng('shuffle')
%
%  E: parameter values at nodes [Nn*1 vector]

% R. Cottereau 10/2008

% constants
if nargin<8; y = []; end
if nargin<9; z = []; end
if nargin<10; dL = []; end
if nargin<11; sd = 'shuffle'; end
if nargin<12; change = []; end
d = 1 + ~isempty(y) + ~isempty(z);
    
% reshaping
if length(L)==1; L = repmat( L, [1 d] ); end
if length(dL)==1; dL = repmat( dL, [1 d] ); end
if size(x,1)==1; x = x'; end
if size(y,1)==1; y = y'; end
if size(z,1)==1; z = z'; end

% deterministic case
if s==0
    k = mu*ones( length(x), max(length(y),1), max(length(z),1), Nmc );
    return
end

% generate unit centered gaussian random field with given correlation 
% structure
if d==1;
    k = makeGauss1D( correl, x, L, dL, Nmc, sd );
elseif d==2;
    k = makeGauss2D( correl, x, y, L, dL, Nmc, sd, change );
elseif d==3;
    k = makeGauss3D( correl, x, y, z, L, dL, Nmc, sd );
end

% use isoprobabilistic transformation to create the random field with
% chosen first-order marginal law, given the gaussian random field above
switch law
    
    % transforming into gamma random field
    case 'gamma'
        k = normcdf( k, 0, 1 );
        k = gaminv( k, mu^2/s, s/mu );
        
    % transforming into uniform random field
    case 'uniform'
        k = normcdf( k, 0, 1 );
        a = sqrt(12*s)/2;
        k = mu + (k-0.5)*(2*a);
        
    % transforming into gauss random field: do nothing
    case {'gauss','normal'}
        k = mu + sqrt(s)*k;
        
    % transforming into lognormal random field
    case 'lognormal'
        c = sqrt(1+s/mu^2);
        k = mu/c*exp(sqrt(log(c^2))*k);
        
    % not implemented yet
    otherwise
        error('this first-order marginal law has not been implemented yet')
        
end

%======= MAKEGAUSS ========
% generator of gaussian random field with unit variance and zero mean
% R. Cottereau 08/2012
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
function gx = makeGauss2D( type, x, y, Lc, dL, Nmc, sd, change )

% constants
xmax = [max(x)-min(x) max(y)-min(y)];
Nfft = 8;

% discretization in fourier space
dk = min( [pi/2./Lc ; pi./xmax] )*.95;
kmax = 4*pi./Lc+dk;
k1 = (0:dk(1):kmax(1));
k2 = (0:dk(2):kmax(2));
Nk = [length(k1) length(k2)];
xk1 = linspace( 0, 2*pi/dk(1), Nfft*Nk(1) );
xk2 = linspace( 0, 2*pi/dk(2), Nfft*Nk(2) );

% construction of power spectral density
switch type
    
    % triangle power spectrum = sinc^2 correlation
    case 'sinc2'
        Sk = tripuls( k1', 4*pi/Lc(1) ) * tripuls( k2, 4*pi/Lc(2) );
        if ~isempty(dL) && any(dL>Lc)
            Sk ( k1>2*pi/dL(1), : ) = 0;
            Sk ( :, k2>2*pi/dL(2) ) = 0;
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
gx = randn( Nk(1), Nk(2), Nmc );
if(~isempty(change) == 1)
    A=change(1);
    B=change(2);
    gx(A:B,A:B, Nmc)= randn( B-A+1, B-A+1, Nmc );
end
%gx(A,B, Nmc)= 100*gx(A,B, Nmc);

% random field in fourier space (isotropic case)
%gk = gx;
gk = fft( fft( gx, [], 1 ), [], 2 );
gk = reshape(reshape(gk,[prod(Nk) Nmc]).*repmat(sqrt(Sk(:)),[1 Nmc] ),[Nk Nmc]);

% return to real space
amp = (Nfft/2/pi/sqrt(2))^2 * prod(Nk.*dk.*Lc);
gx = ifft( ifft( gk, Nfft*Nk(2), 2, 'symmetric' ), ...
                     Nfft*Nk(1), 1, 'symmetric' );
indx = find( xk1>max(x-min(x)), 1 );
indy = find( xk2>max(y-min(y)), 1 );
gx = interp1( xk1(1:indx), gx( 1:indx, 1:indy, : ), x-min(x) );
gx = amp * permute( interp1( xk2(1:indy), ...
                             permute( gx, [2 1 3] ), y-min(y) ), [2 1 3] );

%======= MAKEGAUSS ========
% generator of gaussian random field with unit variance and zero mean
% R. Cottereau 08/2012
function gx = makeGauss3D( type, x, y, z, Lc, dL, Nmc, sd )

% constants
xmax = [max(x)-min(x) max(y)-min(y) max(z)-min(z)];
Nfft = 8;

% discretization in fourier space
dk = min( pi/2./Lc, pi./xmax )*.95;
kmax = 4*pi./Lc+dk;
k1 = (0:dk(1):kmax(1));
k2 = (0:dk(2):kmax(2));
k3 = (0:dk(3):kmax(3));
Nk = [length(k1) length(k2) length(k3)];
xk1 = linspace( 0, 2*pi/dk(1), Nfft*Nk(1) );
xk2 = linspace( 0, 2*pi/dk(2), Nfft*Nk(2) );
xk3 = linspace( 0, 2*pi/dk(3), Nfft*Nk(3) );

% construction of power spectral density
switch type
    
    % triangle power spectrum = sinc^2 correlation
    case 'sinc2'
        Sk = tripuls( k1', 4*pi/Lc(1) ) * tripuls( k2, 4*pi/Lc(2) );
        Sk = reshape( Sk(:) * tripuls( k3, 4*pi/Lc(3) ), Nk );
        if ~isempty(dL) && (dL>Lc)
            Sk ( k1>2*pi/dL, :, : ) = 0;
            Sk ( :, k2>2*pi/dL, : ) = 0;
            Sk ( :, :, k3>2*pi/dL ) = 0;
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
gx = randn( Nk(1), Nk(2), Nk(3), Nmc );

% random field in fourier space (isotropic case)
gk = fft( fft( fft( gx, [], 1 ), [], 2 ), [], 3 );
gk = reshape(reshape(gk,[prod(Nk) Nmc]).*repmat(sqrt(Sk(:)),[1 Nmc] ),[Nk Nmc]);

% return to real space
amp = (Nfft/2/pi/sqrt(2))^3 * prod(Nk.*dk.*Lc);
gx = ifft( ifft( ifft( gk, Nfft*Nk(1), 1, 'symmetric' ), ...
                           Nfft*Nk(2), 2, 'symmetric' ), ...
                           Nfft*Nk(3), 3, 'symmetric' );
indx = find( xk1>xmax(1), 1 );
indy = find( xk2>xmax(2), 1 );
indz = find( xk3>xmax(3), 1 );
gx = interp1( xk1(1:indx), gx( 1:indx, 1:indy, 1:indz, : ), x-min(x) );
gx = permute( interp1( xk2(1:indy), ...
                         permute( gx, [2 1 3 4] ), y-min(y) ), [2 1 3 4] );
gx = amp * permute( interp1( xk3(1:indz), ...
                         permute( gx, [3 2 1 4] ), z-min(z) ), [3 2 1 4] );

