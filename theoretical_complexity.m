function complexity = theoretical_complexity(norm_xStep, norm_L, corrMod, kAdjust, periodMult)


%Call Exemple:
% theoretical_complexity([0.1, 0.1], [100,100], 'gaussian', 10, 1.1)
nDim = size(norm_L,2);
kNTotal = 1;
xNTotal = 1;

if (strcmp(strtrim(corrMod), 'gaussian'))
    corrModChar = 'G';
end

switch corrModChar
    case 'G'
        %ref = pi^nDim;
        %integralExigence = 0.01;
        
        switch nDim
            case 1
                kMax = 6.457; %integralExigence = 0.01
                %                 fun = @(x) exp(-x.^2/(4*pi));
                %                 while(proportion > integralExigence)
                %                     kMax = kMax+0.001;
                %                     trial = quad(fun,0,kMax);
                %                     proportion = 1 - trial/ref;
                %                 end
            case 2
                kMax = 7.035; %integralExigence = 0.01
                %                 fun = @(x,y) exp(-(x.^2+y.^2)/(4*pi));
                %                 while(proportion > integralExigence)
                %                     kMax = kMax+0.001;
                %                     trial = quad2d(fun,0,kMax,0,kMax);
                %                     proportion = 1 - trial/ref;
                %                 end
            case 3
                kMax = 7.355; %integralExigence = 0.01
                %                 fun = @(x,y,z) exp(-(x.^2+y.^2+z.^2)/(4*pi));
                %                 while(proportion > integralExigence)
                %                     kMax = kMax+0.001;
                %                     trial = triplequad(fun,0,kMax,0,kMax,0,kMax);
                %                     proportion = 1 - trial/ref
                %                 end
        end
    otherwise
        error('corMod not implemented')
end

for i = 1:nDim
    kNTotal = kNTotal * (1 + kAdjust*(ceil(periodMult*norm_L(i)*kMax/(2*pi))));
    xNTotal = xNTotal * (1 + floor(norm_L(i)/norm_xStep(i)));
end

complexity = kNTotal * xNTotal * nDim; %When using DGEMM - O(N^3)

end
