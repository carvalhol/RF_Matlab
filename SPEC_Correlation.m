function [C,r]=SPEC_Correlation(parameter1,parameter2,x,y,dir)
% SPEC_CORREL_TRACES to correlate signals
%
% syntax [C,r]=SPEC_Correlation(parameter1,parameter2,x,y,dir)
%
%   parameter1: first parameter
%   parameter2: second parameter
% 
%
%   C: correlation vector following a precised direction
%   r: lag distance

% S. Khazaie 03/2014

% parameter1 and parameter2 should be interpolated before inserting to the 
% fft or ifft 
% Regular grid mesh (over which the interpolation will be done)
xx = min(x):floor(mean(diff(x))):max(x);
yy = min(y):floor(mean(diff(y))):max(y);

[X,Y]=ndgrid(xx,yy);

% valinterp = parameter; % interp2(x,y,parameter,X,Y);
valinterp1 = interp2(x,y,parameter1,X,Y);
valinterp2 = interp2(x,y,parameter2,X,Y);

[n,m] = size(valinterp1);
% r=0:dt:(n-1)*dt;

rho = corrcoef(valinterp1(:),valinterp2(:));
% corcoef = rho(1,2);

if strcmp(dir,'x')
    
    C = ifft( mean(fft(valinterp1-repmat(mean(mean(valinterp1,1),2),[n m]),[],2) ...
        .*conj(fft(valinterp2-repmat(mean(mean(valinterp2,1),2),[n m]),[],2)),1) );
    
    C = C(1:floor(m/2))./C(1);  % (std(valinterp1(:))*std(valinterp2(:))*corcoef)
    
    r = 0 : floor(mean(diff(x))) : floor(mean(diff(x)))*(floor(m/2)-1);
    
elseif strcmp(dir,'y')
    
    C = ifft( mean(fft(valinterp1-repmat(mean(mean(valinterp1,1),2),[n m]),[],1) ...
        .*conj(fft(valinterp2-repmat(mean(mean(valinterp2,1),2),[n m]),[],1)),2) );
    C = C(1:floor(m/2))./C(1);
    r = 0 : floor(mean(diff(y))) : floor(mean(diff(y)))*(floor(m/2)-1);
       
else
    disp('Direction should be either x or y')
end


% if strcmp(dir,'x')
%     
%     C = real(ifft2( mean(fft2(valinterp-repmat(mean(mean(valinterp,1),2),[n m])) ...
%         .*conj(fft2(valinterp-repmat(mean(mean(valinterp,1),2),[n m]))),1) ));
%     C = C(1:floor(n/2))./C(1);
%     r = 0 : floor(mean(diff(x))) : floor(mean(diff(x)))*(floor(n/2)-1);
%     
% elseif strcmp(dir,'y')
%     
%     C = real(ifft2( mean(fft2(valinterp-repmat(mean(mean(valinterp,1),2),[n m])) ...
%         .*conj(fft2(valinterp-repmat(mean(mean(valinterp,1),2),[n m]))),2) ));
%     C = C(1:floor(m/2))./C(1);
%     r = 0 : floor(mean(diff(y))) : floor(mean(diff(y)))*(floor(m/2)-1);
%        
% else
%     disp('Direction should be either x or y')
% end