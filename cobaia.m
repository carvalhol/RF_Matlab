xStep = 0.1;
corrMod = 'G';
kAdjust = 1;
nDim = 3;
nIter = 20;

complexity = zeros(nIter,1);

for i = 1:nIter
    
    normL = zeros(nDim,1);
    normL(:) = i;
    xNTotal = prod(normL)/(xStep^nDim);
    complexity(i) = theoretical_complexity(xNTotal, normL, corrMod, kAdjust);
    
end