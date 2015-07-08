%% Me 1D
% 
% clc
% clear all
% close all
% 
% k     = linspace(0, 50, 1024);
% N     = length(k);
% klong = [k, k(end)+k(2:end-1)];
% 
% A     = exp(-1i*k'.^2);
% Along = exp(-1i*klong'.^2);
% Along(N + 3:end) = 2; %Ignored
% 
% Asym  = [A; conj(flipud(A(2:end-1)))];
% 
% B     = ifft(A,'symmetric');
% Blong = ifft(Along,'symmetric');
% Bsym  = ifft(Asym);
% 
% 
% 
% plot(2*(1:N), B/2, 'b-x',...
%      1:2*N-2, Blong, 'r-x', ...
%      1:2*N-2, Bsym, 'g--')
% 
% figure(1)
% plot(1:2*N-2, abs(Blong), 'b', ...
%      1:2*N-2, abs(Bsym), 'r--')
%  
%  figure(2)
%  plot(1:2*N-2, unwrap(angle(Blong)), 'b', ...
%      1:2*N-2, unwrap(angle(Bsym)), 'r--')
%  
% legendInfo{1} = ['B using symmetric'];
% legendInfo{1} = ['B doubled using symmetric'];
% legendInfo{2} = ['B symmetrized (without symmetric)'];
% legend(legendInfo,'Location','northeast','FontSize',15);

%% Me 2D
% 
% clc
% clear all
% 
% k1     = linspace(0, 50, 1024);
% k2     = linspace(0, 50, 1024);
% k = k1'*k2;
% N     = length(k1);
% k1long = linspace(0, (2*50)-2, (2*1024) -2);
% k2long = linspace(0, (2*50)-2, (2*1024) -2);
% klong = k1long'*k2long;
% 
% A     = exp(-1i*k'.^2);
% Along = exp(-1i*klong'.^2);
% Along(N + 3:end) = 2;
% 
% Asym  = [A, A(:,2:end-1); A(2:end-1,:), A(2:end-1,2:end-1)];
% 
% Asym(N+1:end, :) = conj(flipud(Asym(N+1:end, :)));
% Asym(:, N+1:end) = conj(fliplr(Asym(:, N+1:end)));
% Asym(N+1:end, N+1:end) = conj(Asym(N+1:end, N+1:end));
% 
% B     = ifft2(A,'symmetric');
% Blong = ifft2(Along,'symmetric');
% Bsym  = ifft2(Asym);



% plot(2*(1:N), B/2, 'b-x',...
%      1:2*N-2, Blong, 'r-x', ...
%      1:2*N-2, Bsym, 'g--')
% legendInfo{1} = ['B using symmetric'];
% legendInfo{2} = ['B doubled using symmetric'];
% legendInfo{3} = ['B symmetrized (without symmetric)'];
% legend(legendInfo,'Location','northeast','FontSize',15);

% %Modified Vector
% Aplus = [A; conj(flipdim(A,1))]
% 
% %iffts
% B = ifft( A, [], 1, 'symmetric' )
% Bplus = ifft( Aplus, [], 1)

%% Regis 1D

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

%% Regis 2D

clc
clear all
close all

%Original matrix
Ainit = [zeros(685,1); sin((0:1:628)'/100); zeros(686,1)];
N = length(Ainit);
Ainit = Ainit*Ainit';
A = fft2(Ainit, N, N);

%Modified Vector
Nhalf = floor(N/2)+1;
Ahalf = A(1:Nhalf, 1:Nhalf);

%Copy
Aplus = [Ahalf, Ahalf(:,2:end-1); ...
         Ahalf(2:end-1,:), Ahalf(2:end-1,2:end-1)];
     
%Hermitian Conjugate
Aplus(Nhalf+1:end, :) = -conj(flipud(Aplus(Nhalf+1:end, :)));
Aplus(:, Nhalf+1:end) = -conj(fliplr(Aplus(:, Nhalf+1:end)));
 
%iFFTs
%B = ifft(A);
Bsym = ifft2( A, 'symmetric' );
Bplus = ifft2( Aplus);

figure; surf(Bsym); shading flat; view(2); colorbar; title('B sym')
figure; surf(real(Bplus)); shading flat; view(2); colorbar; title('B+')
% 
% %Original matrix
% N = 100;
% Ainit = [zeros(N-1,1); sin((0:N)'*2*pi/N); zeros(N,1)];
% N = length(Ainit);
% Ainit = Ainit*Ainit';
% figure; surf(Ainit); shading flat; view(2); colorbar; title('A')
% A = fft2(Ainit);
% 
% %Modified Vector
% Nhalf = floor(N/2)+1;
% 
% %Copy
% Aplus = A;
% Aplus( Nhalf+1:end, 2:Nhalf ) = -conj(flipud(A( 2:Nhalf-1, 2:Nhalf )));
% Aplus( 2:Nhalf, Nhalf+1:end ) = -conj(fliplr(A( 2:Nhalf, 2:Nhalf-1 )));
% Aplus( Nhalf+1:end, Nhalf+1:end ) = rot90(A( 2:Nhalf-1, 2:Nhalf-1 ),2);
% 
% %iFFTs
% Bsym = ifft2( A, 'symmetric' );
% figure; surf(Bsym); shading flat; view(2); colorbar; title('B sym')
% Bplus = ifft2(Aplus);
% figure; surf(real(Bplus)); shading flat; view(2); colorbar; title('B+')



