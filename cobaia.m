xStepBase = 0.1;
corrMod = 'gaussian';
kAdjust = 1;
nDim = 2;
nIter = 4;
periodMult = 1.1;

complexity = zeros(nIter,1);
kNTotal = zeros(nIter,1);
xNTotal = zeros(nIter,1);
ratio = zeros(nIter,1);

normL = zeros(nDim,1);
xStep = zeros(nDim,1);

normL(:) = 100;
xStep(:) = xStepBase;
xAreaTotal = 10^(6);

for i = 1:nIter
    
    normL(1) = 10^i;
    normL(2) = round(xAreaTotal/normL(1))    

    ratio(i) = normL(2)/normL(1)
    [complexity(i),xNTotal(i),kNTotal(i)] = theoretical_complexity(xStep, normL, corrMod, kAdjust, periodMult)
    
end
