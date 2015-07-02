%% Me

clc
clear all

k     = linspace(0, 50, 1024);
N     = length(k);
klong = linspace(0, (2*50)-2, (2*1024) -2);

A     = exp(-1i*k'.^2);
Along = exp(-1i*klong'.^2);
%Along(N + 3:end) = 2;

Asym  = [A; conj(flipud(A(2:end-1)))];

B     = ifft(A,'symmetric');
Blong = ifft(Along,'symmetric');
Bsym  = ifft(Asym);



plot(2*(1:N), B/2, 'b-x',...
     1:2*N-2, Blong, 'r-x', ...
     1:2*N-2, Bsym, 'g--')
legendInfo{1} = ['B using symmetric'];
legendInfo{2} = ['B doubled using symmetric'];
legendInfo{3} = ['B symmetrized (without symmetric)'];
legend(legendInfo,'Location','northeast','FontSize',15);

% %Modified Vector
% Aplus = [A; conj(flipdim(A,1))]
% 
% %iffts
% B = ifft( A, [], 1, 'symmetric' )
% Bplus = ifft( Aplus, [], 1)

%% Regis
%Original vector
% Ainit = [zeros(100,1); sin((0:0.1:2*pi)'); zeros(100,1)];
% A = fft( Ainit );
% N = length(A);
% %Modified Vector
% Aplus = A;
% Aplus(end-floor(N/2)+1:end) = conj(fliplr(A(1:floor(N/2))));
% figure; plot(1:N,abs(A),'b-x')
% figure; plot(1:N,abs(Aplus),'r-x');
% 
% %iffts
% B = ifft( A, [], 1, 'symmetric' );
% Bplus = ifft( Aplus(1:end-1), [], 1);
% C = ifft( A, [], 1 );
% figure; plot(1:N,B,'b-x',1:N-1,real(Bplus),'r-x')
%figure; plot(1:N,B,'b-x')
%figure; plot(1:N-1,real(Bplus),'r-x');



% Aconj = conj(flipdim(A,2))
% Aplus = [A(:,1:2),conj(flipdim(A(:,1:2),2))]
% 
% B = ifft( A, [], 2, 'symmetric' )
% Bplus = ifft( Aplus, [], 2)

% NFFT = 128;
% x = randn(NFFT,1);
% H = zeros(NFFT,1);
% H(10:20) = 1;
% y2 = ifft(H.*fft(x), 'symmetric');
% 
% % H(end-20+2:end-10+2) = 1;    % Other half
% y  = ifft(H.*fft(x));
% 
% y2-y

