function [C, lc]=corr_RF(RF_1,RF_2, dir, xStep_dir)

% parameter1 and parameter2 should be interpolated before inserting to the 
% fft or ifft 
% Regular grid mesh (over which the interpolation will be done)
%xx = min(x):floor(mean(diff(x))):max(x);
%yy = min(y):floor(mean(diff(y))):max(y);

%[x_Uni,y_Uni]=ndgrid(xx,yy);

%RF_Uni_1 = interp2(x,y,RF1,x_Uni,y_Uni);
%RF_Uni_2 = interp2(x,y,RF2,x_Uni,y_Uni);

mean_1 = mean(RF_1(:));
mean_2 = mean(RF_2(:));

nDim = 2;
dims = 1:nDim;
meanDims = dims(setdiff(1:end,dir));

C = fft(RF_1-mean_1,[],dir).*conj(fft(RF_2-mean_2,[],dir));
%C = C./numel(C);
   
for i = 1 : size(meanDims)
    C = mean(C, meanDims(i));
end

C = ifft(C);
C = C/C(1);
lc = trapz(C)*xStep_dir*(numel(C)-1)/(numel(C));
%C = C.*numel(C);
C = C(1:floor(numel(C)/2))./C(1);

if(size(C,1) > size(C,2))
    C = C';
end
%C = C(1:floor(size(C)/2));
    
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