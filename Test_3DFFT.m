function Y = Test_3DFFT(x,y,z,corrModx,corrMody,corrModz, Lx, Ly, Lz, delta)

%In space domain
hx = x(2) - x(1);
hy = y(2) - y(1);
hz = z(2) - z(1);
x  = min(x):hx:max(x);
y  = min(y):hy:max(y);
z  = min(z):hz:max(z);
Nx = 2*numel(x) ;
Ny = 2*numel(y) ;
Nz = 2*numel(z) ;

%In wave-number domain
Tx = max(x) - min(x) ; dwx = 2*pi/Tx ; wx = (0:Nx-1)*dwx ; %Nx or Nx/2 ??
Ty = max(y) - min(y) ; dwy = 2*pi/Ty ; wy = (0:Ny-1)*dwy ;
Tz = max(z) - min(z) ; dwz = 2*pi/Tz ; wz = (0:Nz-1)*dwz ;

%Random Variables
phik   = rand(Nx,Ny,Nz);
gammak = randn(Nx,Ny,Nz) ;

%Spectrum
Sk  = zeros(Nx,Ny,Nz);
Skx = dwx.*correlationModel(wx, corrModx , Lx, delta) ;
Sky = dwy.*correlationModel(wy, corrMody , Ly, delta) ;
Skz = dwz.*correlationModel(wz, corrModz , Lz, delta) ;

Sk = Lx*Ly*Lz * Skx' * Sky; 

Sk = repmat( Sk, [1 1 Nz]) .* ... 
     repmat( reshape(Skz, ... 
     [1 1 Nz]), [Nx Ny 1]);
 
 %Phase and amplitudes
Dk     = gammak.*sqrt(Sk).*exp(2*pi*phik);

 
 %START New Spectrum Organization
 L    = [Lx, Ly, Lz]';
 xMax = [Tx, Ty, Tz]';
 nDim = 3;
 kMax = zeros(nDim, 1);
 kMin = zeros(nDim, 1);
 %kMax(:) = 7.355;

 kStepMax = 2*pi./(xMax);
 kMax  = kStepMax.*[(Nx-1), (Ny-1), (Nz-1)]';
 kStep = kMax./(ceil(kMax./kStepMax)+1);
 kNStep = kMax./kStep;
 %kVec  = 0:kStep:kMax;
 kPoints = zeros(prod(kNStep),nDim);
 
 for pos = 1:size(kPoints,1)
     kPoints(pos,:) = get_perm(pos, kMax, kNStep', kMin);
 end
 
 
 
 Sk2 = prod(kStepMax).*correlationModel2(kPoints, corrModx , L, delta) ;
 
 phik2   = rand(Nx*Ny*Nz,1);
 gammak2 = randn(Nx*Ny*Nz,1);
 
 Dk2 = gammak2.*sqrt(Sk2).*exp(complex(0,1)*2*pi*phik2);
 
 %END New Spectrum Organization

% %Field computation

Y = ifft( ifft( ifft( Dk, [], 3, 'symmetric' ), ... 
    [], 2, 'symmetric' ), ... 
    [], 1, 'symmetric' ); 
Y = Y(1:Nx/2, 1:Ny/2, 1:Nz/2) ...
    * sqrt(2*pi/Lx)*Nx ... 
    * sqrt(2*pi/Ly)*Ny ... 
    * sqrt(2*pi/Lz)*Nz;

% Y = real(numel(DkPlus)*ifftn(DkPlus)); %Would be different if we use Y = numel(DkPlus)*ifft2(DkPlus,'symmetric');
% Y = Y(Nx:2*Nx, Ny:2*Ny, Nz:2*Nz); %Why this quadrant?
% Y = Y(1:2:end , 1:2:end, 1:2:end); %Is this only cause you put h/2 in the begining?
  
end

function Sk = correlationModel( wk, model,L ,delta )
switch lower(model)
    case 'cardinalsine2'
        eta = pi*L ;
        Sk  = 1/(eta*pi)*tripuls( wk,2/eta);
        
    case 'gaussian'
        %eta = L/sqrt(pi);
        %Sk = 1/(2*pi)*exp(-eta^2/4*wk.^2);
        Sk = 1/(2*pi)*sqrt(pi/log(1/delta))*L*exp( -wk.^2*L^2/(4*log(1/delta)) );
        
    case 'exponential'
        Sk = 1/pi*(L/2)./(1+(L/2*wk).^2);
        
    otherwise
        error('this correlation model does not exist')
        
end
end

function SkVec = correlationModel2( wk, model,L ,delta )
switch lower(model)

    case 'gaussian'
        %eta = L/sqrt(pi);
        %Sk = 1/(2*pi)*exp(-eta^2/4*wk.^2);
        SkVec = ones(size(wk,1),1);
        for i = 1:3
            SkVec = SkVec.*(1/(2*pi)*sqrt(pi/log(1/delta))*L(i)*exp( -wk(:,i).^2*L(i)^2/(4*log(1/delta)) ));
        end
end
end

function M_out = flip(M_in, boolVector)

M_out = M_in;

for i = 1:numel(boolVector)
    if(boolVector(i))
        M_out = flipdim(M_out, i);
    end
end

end