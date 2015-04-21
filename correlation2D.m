function corr = correlation2D(H , x , y ) 
S = size(H) ;
Nx = S(1)  ;
Ny = S(2) ;
Nmc = S(3) ;

 

 
% Calcul des moyennes 
% m0 = mean(H(x,y,:)) ;
m = zeros(Nx,Ny) ;
for k = 1:Nmc
    m = m + 1/Nmc*H(:,:,k) ;
end

 
% Calcul des écarts-types
sigma0 = std(H(x,y,:)) ;
sigma = zeros(Nx,Ny) ; 
for k = 1:Nmc
    sigma = sigma + 1/Nmc*(H(:,:,k)-m).^2 ; % Contribution pour la variance
end
sigma = sqrt(sigma) ;

 
%Calcul de la corrélation au point (x,y) 
corr = zeros(Nx,Ny) ; 
for k = 1:Nmc 
    corr = corr + 1/Nmc*H(x,y,k)*H(:,:,k) ;
end

 
corr = corr/sigma0./sigma;