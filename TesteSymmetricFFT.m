%% 1D

% clc
% clear all
% close all

% %Original vector
% Ainit = [zeros(685,1); sin((0:1:628)'/100); zeros(686,1)];
% A = fft( Ainit );
% N = length(A);
% %Modified Vector
% Aplus = A;
% Aplus(end-floor(N/2)+1:end) = conj(flipud(A(1:floor(N/2))));
% figure; plot(1:N,abs(A),'b-x',1:N,abs(Aplus),'r-x');
% 
% %iffts
% B = ifft(A);
% B = ifft( A, [], 1, 'symmetric' );
% Bplus = ifft( Aplus(1:end-1), [], 1);
% figure; plot(1:N,B,'b-x',1:N-1,real(Bplus),'r-x');

%% 2D

% clc
% clear all
% close all
% 
% %Original matrix
% Ainit = [zeros(685,1); sin((0:1:628)'/100); zeros(686,1)];
% N = length(Ainit);
% Ainit = Ainit*Ainit';
% A = fft2(Ainit, N, N);
% 
% %Modified Vector
% Nhalf = floor(N/2)+1;
% Ahalf = A(1:Nhalf, 1:Nhalf);
% 
% %Copy
% Aplus = [Ahalf, Ahalf(:,2:end-1); ...
%          Ahalf(2:end-1,:), Ahalf(2:end-1,2:end-1)];
%      
% %Hermitian Conjugate
% Aplus(Nhalf+1:end, :) = -conj(flipud(Aplus(Nhalf+1:end, :)));
% Aplus(:, Nhalf+1:end) = -conj(fliplr(Aplus(:, Nhalf+1:end)));
%  
% %iFFTs
% %B = ifft(A);
% Bsym = ifft2( A, 'symmetric' );
% Bplus = ifft2( Aplus);
% 
% figure; surf(Bsym); shading flat; view(2); colorbar; title('B sym')
% figure; surf(real(Bplus)); shading flat; view(2); colorbar; title('B+')

%% 3D

clc
clear all
close all

%Original matrix
%A1D = [zeros(5,1); sin((0:1:628)'/100); zeros(6,1)];
A1D = [1:1:100]'; %A1 should be even
N = length(A1D);
Ainit = A1D*A1D';

Ainit = repmat( Ainit, [1 1 N]) .* ... 
        repmat( reshape(A1D, ... 
        [1 1 N]), [N N 1]);
 
A = fftn(Ainit);

%Modified Vector
Nhalf = floor(N/2)+1;
Ahalf = A(1:Nhalf, 1:Nhalf, 1:Nhalf);

%Copy
Aplus = cat(1,...
        cat(2, Ahalf, Ahalf(:,2:end-1,:)), ...
        cat(2, Ahalf(2:end-1,:,:), Ahalf(2:end-1,2:end-1,:))...
        );
         
Aplus = cat(3, Aplus, Aplus(:,:,2:end-1));

% for i=1:N
%     Aplus(end+1-i,i,i)
% end
      
 %Hermitian Conjugate
 Aplus(Nhalf+1:end, :,:) = -conj(flipdim(Aplus(Nhalf+1:end, :,:),1));
 Aplus(:, Nhalf+1:end,:) = -conj(flipdim(Aplus(:, Nhalf+1:end,:),2));
 Aplus(:,:, Nhalf+1:end) = -conj(flipdim(Aplus(:,:, Nhalf+1:end),3));


%iFFTs
%%B = ifft(A);
Bsym = ifftn( A, 'symmetric' );
Bplus = ifftn( Aplus);

[x,y,z] = meshgrid(1:N,1:N,1:N);
scatter3(x(:),y(:),z(:),5,real(Aplus(:)))

% figure; surf(Bsym); shading flat; view(2); colorbar; title('B sym')
% figure; surf(real(Bplus)); shading flat; view(2); colorbar; title('B+')