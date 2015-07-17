function Y = ranfomfieldsgeneration3DFFT(x,y,z,corrModx,corrMody,corrModz, Lx, Ly, Lz, delta)

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

% zk = sqrt( -log( rand( Nx, Ny, Nz, Nmc ) ) );
% phik = 2*pi * rand( Nx, Ny, Nz, Nmc );
% 
% Sk = repmat( sqrt( Sk ), [1 1 1 Nmc] );
% gk = Sk .* zk .* exp( complex(0,1)*phik );

% g = ifft( ifft( ifft( gk, [], 3, 'symmetric' ), ... 
%     [], 2, 'symmetric' ), ... 
%     [], 1, 'symmetric' ); 
% g = g(1:Nx/2,1:Ny/2,1:Nz/2,:) * sqrt(2*pi/L(1))*Nx ... 
%     * sqrt(2*pi/L(2))*Ny ... 
%     * sqrt(2*pi/L(3))*Nz;

% for i = 1:Nx
%     for j = 1:Ny
%         for k = 1:Nz
%             Sk(i,j,k) = Skx(i)*Sky(j)*Skz(k);
%         end
%     end
% end

%Phase and amplitudes
Dk     = gammak.*sqrt(Sk).*exp(2*pi*phik);
% DkPlus = cat(3, ...
%     cat(2, ...
%     cat(1,flip(Dk,[1,1,1]), flip(Dk,[0,1,1])), ...
%     cat(1,flip(Dk,[1,0,1]), flip(Dk,[0,0,1])) ...
%     ), ...
%     cat(2, ...
%     cat(1,flip(Dk,[1,1,0]), flip(Dk,[0,1,0])), ...
%     cat(1,flip(Dk,[1,0,0]), flip(Dk,[0,0,0])) ...
%     )  ...
%     );

%Only To verify symetry
% shift = 2;
% 
% a = [1-shift,1-shift,1-shift];
% DkPlus(Nx+a(1), Ny+a(2), Nz+a(3))
% a = [1-shift,1-shift,shift];
% DkPlus(Nx+a(1), Ny+a(2), Nz+a(3))
% a = [1-shift,shift,1-shift];
% DkPlus(Nx+a(1), Ny+a(2), Nz+a(3))
% a = [1-shift,shift,shift];
% DkPlus(Nx+a(1), Ny+a(2), Nz+a(3))
% a = [shift,1-shift,1-shift];
% DkPlus(Nx+a(1), Ny+a(2), Nz+a(3))
% a = [shift,1-shift,shift];
% DkPlus(Nx+a(1), Ny+a(2), Nz+a(3))
% a = [shift,shift,1-shift];
% DkPlus(Nx+a(1), Ny+a(2), Nz+a(3))
% a = [shift,shift,shift];
% DkPlus(Nx+a(1), Ny+a(2), Nz+a(3))


% %Field computation

% Y = ifft( ifft( ifft( Dk, [], 3, 'symmetric' ), ... 
%     [], 2, 'symmetric' ), ... 
%     [], 1, 'symmetric' ); 

N = [Nx, Ny, Nz];

Nhalf = floor(N/2)+1;
Dkhalf = Dk(1:Nhalf(1), 1:Nhalf(2), 1:Nhalf(3));

%Copy
Dkplus = cat(1,...
        cat(2, Dkhalf, Dkhalf(:,2:end-1,:)), ...
        cat(2, Dkhalf(2:end-1,:,:), Dkhalf(2:end-1,2:end-1,:))...
        );
         
Dkplus = cat(3, Dkplus, Dkplus(:,:,2:end-1));

%Hermitian Conjugate
Dkplus(Nhalf+1:end, :,:) = -conj(flipdim(Dkplus(Nhalf+1:end, :,:),1));
Dkplus(:, Nhalf+1:end,:) = -conj(flipdim(Dkplus(:, Nhalf+1:end,:),2));
Dkplus(:,:, Nhalf+1:end) = -conj(flipdim(Dkplus(:,:, Nhalf+1:end),3));

%Field Generation

% Y = ifft( ifft( ifft( Dk, [], 3, 'symmetric' ), ... 
%     [], 2, 'symmetric' ), ... 
%     [], 1, 'symmetric' ); 

%Y = ifftn(Dk, [],'symmetric' ); 


%START Test Without Symmetric
Y = real(ifftn(Dk)); 
%END Test Without Symmetric

%Normalization
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

function M_out = flip(M_in, boolVector)

M_out = M_in;

for i = 1:numel(boolVector)
    if(boolVector(i))
        M_out = flipdim(M_out, i);
    end
end

end