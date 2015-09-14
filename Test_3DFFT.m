function Y = Test_3DFFT(x,y,z,corrModx,corrMody,corrModz, Lx, Ly, Lz, delta)

%In space domain
hx = x(2) - x(1);
hy = y(2) - y(1);
hz = z(2) - z(1);
x  = min(x):hx:max(x);
y  = min(y):hy:max(y);
z  = min(z):hz:max(z);
% Nx = 2*numel(x) ;
% Ny = 2*numel(y) ;
% Nz = 2*numel(z) ;
Nx = numel(x) ;
Ny = numel(y) ;
Nz = numel(z) ;


adjust = 1.1;

%In wave-number domain
Tx = max(x) - min(x) ; dwx = 2*pi/Tx/adjust ; wx = (0:Nx-1)*dwx ; %Nx or Nx/2 ??
Ty = max(y) - min(y) ; dwy = 2*pi/Ty/adjust ; wy = (0:Ny-1)*dwy ;
Tz = max(z) - min(z) ; dwz = 2*pi/Tz/adjust ; wz = (0:Nz-1)*dwz ;

%Random Variables
phik   = rand(Nx,Ny,Nz);
gammak = randn(Nx,Ny,Nz) ;

%Spectrum
Sk  = zeros(Nx,Ny,Nz);
Skx = correlationModel(wx, corrModx , Lx, delta) ;
Sky = correlationModel(wy, corrMody , Ly, delta) ;
Skz = correlationModel(wz, corrModz , Lz, delta) ;

%Sk = Lx*Ly*Lz * dwx*dwy*dwz * Skx' * Sky; 
Sk = Skx' * Sky; 

Sk = repmat( Sk, [1 1 Nz]) .* ... 
     repmat( reshape(Skz, ... 
     [1 1 Nz]), [Nx Ny 1]);
 
 fprintf('Sk---------\n');
 
 i = 1; Sk(i,i,i)
 i = 2; Sk(i,i,i)
 i = 3; Sk(i,i,i)
 i = 4; Sk(i,i,i)
 i = 5; Sk(i,i,i)

 
 %Phase and amplitudes
%Dk     = gammak.*sqrt(Sk).*exp(2*pi*phik);
Dk     = sqrt(Sk);

% ContructingSym
DkSym = zeros(2*size(Dk)-2);

%Symmetrization
N = [Nx; Ny; Nz];

kS  = ones(3,1)*2;
kE  = N-1;
kSc = N+ 1;
kEc = 2*N-2;
        
%Copy
DkSym(1:N(1), 1:N(2), 1:N(3)) = Dk(:,:,:);
 
 DkSym(kSc(1):kEc(1)  , 1:N(2), 1:N(3)) = Dk(kS(1):kE(1),:,:);
 DkSym(1:N(1), kSc(2):kEc(2)  , 1:N(3)) = Dk(:,kS(2):kE(2),:);
 DkSym(1:N(1), 1:N(2), kSc(3):kEc(3)  ) = Dk(:,:,kS(3):kE(3));
 
 DkSym(kSc(1):kEc(1)  , kSc(2):kEc(2)  , 1:N(3)) = Dk(kS(1):kE(1), kS(2):kE(2), :          ); %[-1,-1, 1]
 DkSym(1:N(1), kSc(2):kEc(2)  , kSc(3):kEc(3)  ) = Dk(:          , kS(2):kE(2), kS(3):kE(3)); %[ 1,-1,-1]
 DkSym(kSc(1):kEc(1)  , 1:N(2), kSc(3):kEc(3)  ) = Dk(kS(1):kE(1), :          , kS(3):kE(3)); %[-1, 1,-1]
 
 DkSym(kSc(1):kEc(1)  , kSc(2):kEc(2)  , kSc(3):kEc(3)  ) = Dk(kS(1):kE(1), kS(2):kE(2), kS(3):kE(3)); %[-1,-1,-1]
 
 %Hermitian Conjugate
 DkSym(kSc(1):kEc(1), :, :) = -(DkSym(kEc(1):-1:kSc(1), :, :));
 DkSym(:, kSc(2):kEc(2), :) = -(DkSym(:, kEc(2):-1:kSc(2), :));
 DkSym(:, :, kSc(3):kEc(3)) = -(DkSym(:, :, kEc(3):-1:kSc(3)));

 fprintf('Hermitian Conjugate\n');
 
 i = 1; DkSym(i,i,i)
 i = 2; DkSym(i,i,i)
 i = 3; DkSym(i,i,i)
 i = 4; DkSym(i,i,i)
 i = 5; DkSym(i,i,i)
 
 fprintf('Hermitian Conjugate OTHER\n');
 i = 40; DkSym(i,i,i)
 i = 39; DkSym(i,i,i)
 i = 38; DkSym(i,i,i)
 i = 37; DkSym(i,i,i)
 i = 36; DkSym(i,i,i)
 
% %Field computation

% Y = ifft( ifft( ifft( DkSym, [], 3, 'symmetric' ), ... 
%      [], 2, 'symmetric' ), ... 
%      [], 1, 'symmetric' )*Nx*Ny*Nz; 
 
%  Y = ifft( ifft( ifft( DkSym, [], 3 ), ... 
%      [], 2 ), ... 
%      [], 1 ); 
 Y = ifft( ifft( ifft( DkSym, [], 3, 'symmetric' ), ... 
                [], 2, 'symmetric' ), ... 
                [], 1, 'symmetric' ); 

 
fprintf('Generated Field\n'); 
 i = 1; Y(i,i,i)
 i = 2; Y(i,i,i)
 i = 3; Y(i,i,i)
 i = 4; Y(i,i,i)
 i = 5; Y(i,i,i)
 
 fprintf('Generated Field OTHER \n'); 
 i = 40; Y(i,i,i)
 i = 39; Y(i,i,i)
 i = 38; Y(i,i,i)
 i = 37; Y(i,i,i)
            
 Y = Y(1:Nx, 1:Ny, 1:Nz);



 
% Y = Y(1:Nx, 1:Ny, 1:Nz) ...
%     * sqrt(2*pi/Lx) ... 
%     * sqrt(2*pi/Ly) ... 
%     * sqrt(2*pi/Lz);


% 
% Y = ifft( ifft( ifft( Dk, [], 3, 'symmetric' ), ... 
%      [], 2, 'symmetric' ), ... 
%      [], 1, 'symmetric' ); 
% Y = Y(1:Nx/2, 1:Ny/2, 1:Nz/2) ...
%     * sqrt(2*pi/Lx)*Nx ... 
%     * sqrt(2*pi/Ly)*Ny ... 
%     * sqrt(2*pi/Lz)*Nz;


% Y = ifft( ifft( ifft( DkSym, [], 3 ), ... 
%     [], 2), ... 
%     [], 1); 
% Y = Y(1:Nx/2, 1:Ny/2, 1:Nz/2) ...
%     * sqrt(2*pi/Lx)*Nx ... 
%     * sqrt(2*pi/Ly)*Ny ... 
%     * sqrt(2*pi/Lz)*Nz;

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
        %OLD Sk = 1/(2*pi)*sqrt(pi/log(1/delta))*L*exp( -wk.^2*L^2/(4*log(1/delta)) );
        Sk = exp(-wk.^2/(4.0d0*pi));
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