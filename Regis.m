%% 1D
clc
clear all
close all

%Original matrix
N = 100;
Ainit = [zeros(N-1,1); sin((0:N)'*2*pi/N); zeros(N,1)];
N = length(Ainit);
Ainit = Ainit*Ainit';
figure; surf(Ainit); shading flat; view(2); colorbar; title('Ainit')
A = fft2(Ainit);
figure; surf(real(A)); shading flat; view(2); colorbar; title('A')

%Modified Vector
Nhalf = floor(N/2)+1;
%A(Nhalf:end,:) = 0.0;
%A(:, Nhalf:end) = 0.0;

%Copy
Aplus = A;
Aplus( Nhalf+1:end, 2:Nhalf ) = -conj(flipud(A( 2:Nhalf-1, 2:Nhalf )));
Aplus( 2:Nhalf, Nhalf+1:end ) = -conj(fliplr(A( 2:Nhalf, 2:Nhalf-1 )));
Aplus( Nhalf+1:end, Nhalf+1:end ) = rot90(A( 2:Nhalf-1, 2:Nhalf-1 ),2);
figure; surf(real(Aplus)); shading flat; view(2); colorbar; title('Aplus')

%iFFTs
Bsym = ifft2( A, 'symmetric' );
!Bsym = ifft2( A);
figure; surf(Bsym); shading flat; view(2); colorbar; title('B sym')
Bplus = ifft2(Aplus);
figure; surf(real(Bplus)); shading flat; view(2); colorbar; title('B+')