function Y = ranfomfieldsgeneration2DFFT(x,y,corrModx,corrMody , Lx, Ly, delta) 

%In space domain
hx = x(2) - x(1); 
hy = y(2) - y(1);
x  = min(x):hx/2:max(x); %Is hx/2 mandatory? Only for aliaising?
y  = min(y):hy/2:max(y);
Nx = size(x) ; Nx = Nx(2);
Ny = size(y) ; Ny = Ny(2);

%In wave-number domain
Tx = max(x) - min(x ) ; dwx = 2*pi/Tx ; wx = (0:Nx-1)*dwx ;
Ty = max(y) - min(y) ; dwy = 2*pi/Ty ; wy = (0:Ny-1)*dwy ;

%Random Variables
phik   = rand(Nx,Ny); 
gammak = randn(Nx,Ny) ;

%Spectrum
Skx = dwx.*correlationModel(wx, corrModx , Lx, delta) ;
Sky = dwy.*correlationModel(wy, corrMody , Ly, delta) ;
Sk  = Skx'*Sky ;

%Phase and amplitudes
Dk     = 2*gammak.*sqrt(Sk).*exp(complex(0,1)*2*pi*phik);
DkPlus = [flipdim(Dk,2)           , Dk; ...
          flipdim(flipdim(Dk,1),2), flipdim(Dk,1)]; %What is the influence of the order in Dkplus construction 

%Field computation
Y = real(numel(DkPlus)*ifft2(DkPlus)); %Would be different if we use Y = numel(DkPlus)*ifft2(DkPlus,'symmetric');
Y = Y(Nx:2*Nx,1:Ny+1); %Why this quadrant?
Y = Y(1:2:end , 1:2:end); %Is this only cause you put h/2 in the begining? 

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