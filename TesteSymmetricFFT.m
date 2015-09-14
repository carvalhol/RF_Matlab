
%% 1D
% 
% clc
% clear all
% close all
% 
% 
% 
% %% First Files
% fid = fopen('FFT1D/beforeFFT');
% befFFT = fscanf(fid,'%f %f',[2, Inf]);
% befFFT = befFFT';
% fclose(fid);
% bef = befFFT(:,1);
% 
% fid = fopen('FFT1D/afterFFT');
% aftFFT = fscanf(fid,'%f %f',[2, Inf]);
% aftFFT = aftFFT';
% fclose(fid);
% aft = aftFFT(:,1) + 1i * aftFFT(:,2);
% 
% fid = fopen('FFT1D/afteriFFT');
% aftiFFT = fscanf(fid,'%f %f',[2, Inf]);
% aftiFFT = aftiFFT';
% fclose(fid);
% aftI = aftiFFT(:,1);
% 
% %% Second Files
% 
% fid = fopen('FFT1D/beforeFFT_2');
% befFFT_2 = fscanf(fid,'%f %f',[1, Inf]);
% befFFT_2 = befFFT_2';
% fclose(fid);
% bef = befFFT_2(:,1);
% 
% fid = fopen('FFT1D/afterFFT_2');
% aftFFT_2 = fscanf(fid,'%f %f',[2, Inf]);
% aftFFT_2 = aftFFT_2';
% fclose(fid);
% aft_2 = aftFFT_2(:,1) + 1i * aftFFT_2(:,2);
% 
% fid = fopen('FFT1D/afteriFFT_2');
% aftiFFT_2 = fscanf(fid,'%f %f',[1, Inf]);
% aftiFFT_2 = aftiFFT_2';
% fclose(fid);
% aftI_2 = [aftiFFT_2(:,1)];
% 
% %% Original vector
% N = 2000;
% N_h = N/2;
% 
% Ainit = [zeros(685,1); sin((0:1:628)'/100); zeros(686,1)];
% A = fft( Ainit );
% A_h = A(1:N_h);
% 
% %Comparisons
% figure; plot(1:N,Ainit(1:N),'b-x',1:N,befFFT_2(1:N,1),'r-x');
% figure; plot(1:N,abs(A(1:N)),'b-x',1:N,abs(aft_2(1:N)),'r-x');
% figure; plot(1:N,Ainit(1:N),'b-x',1:N,real(aftiFFT_2(1:N))./N,'r-x');

%% 2D

clc
clear all
close all

%Original matrix
Ainit = [zeros(685,1); sin((0:1:628)'/100); zeros(686,1)];
N = length(Ainit);
Ainit = Ainit*Ainit';
figure; surf(Ainit); shading flat; view(2); colorbar; title('A init')

A = fft2(Ainit, N, N);
figure; surf(abs(A)); shading flat; view(2); colorbar; title('A')



% A(end/2+1:end,:) = [];
% A(:,end/2+1:end) = [];

% A(1:end/2,:) = [];
% A(:,1:end/2) = [];
%figure; surf(abs(A)); shading flat; view(2); colorbar; title('A-')

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
 
%iFFTs
%B = ifft(A);
Bsym = ifft2( A, 'symmetric' );
% Bplus = ifft2( Aplus);

figure; surf(Bsym); shading flat; view(2); colorbar; title('B sym')
%figure; surf(real(Bplus)); shading flat; view(2); colorbar; title('B+')

%% 3D
% 
% clc
% clear all
% close all
% 
% Original matrix
% A1D = [zeros(5,1); sin((0:1:628)'/100); zeros(6,1)];
% A1D = [1:1:100]'; %A1 should be even
% N = length(A1D);
% Ainit = A1D*A1D';
% 
% Ainit = repmat( Ainit, [1 1 N]) .* ... 
%         repmat( reshape(A1D, ... 
%         [1 1 N]), [N N 1]);
%  
% A = fftn(Ainit);
% 
% Modified Vector
% Nhalf = floor(N/2)+1;
% Ahalf = A(1:Nhalf, 1:Nhalf, 1:Nhalf);
% 
% Copy
% Aplus = cat(1,...
%         cat(2, Ahalf, Ahalf(:,2:end-1,:)), ...
%         cat(2, Ahalf(2:end-1,:,:), Ahalf(2:end-1,2:end-1,:))...
%         );
%          
% Aplus = cat(3, Aplus, Aplus(:,:,2:end-1));
% 
% for i=1:N
%     Aplus(end+1-i,i,i)
% end
%       
%  Hermitian Conjugate
%  Aplus(Nhalf+1:end, :,:) = -conj(flipdim(Aplus(Nhalf+1:end, :,:),1));
%  Aplus(:, Nhalf+1:end,:) = -conj(flipdim(Aplus(:, Nhalf+1:end,:),2));
%  Aplus(:,:, Nhalf+1:end) = -conj(flipdim(Aplus(:,:, Nhalf+1:end),3));
% 
% 
% iFFTs
% %B = ifft(A);
% Bsym = ifftn( A, 'symmetric' );
% Bplus = ifftn( Aplus);
% 
% [x,y,z] = meshgrid(1:N,1:N,1:N);
% scatter3(x(:),y(:),z(:),5,real(Aplus(:)))
% 
% figure; surf(Bsym); shading flat; view(2); colorbar; title('B sym')
% figure; surf(real(Bplus)); shading flat; view(2); colorbar; title('B+')