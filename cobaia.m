% xStep = 0.1;
% corrMod = 'G';
% kAdjust = 1;
% nDim = 3;
% nIter = 20;
% 
% complexity = zeros(nIter,1);
% 
% for i = 1:nIter
%     
%     normL = zeros(nDim,1);
%     normL(:) = i;
%     xNTotal = prod(normL)/(xStep^nDim);
%     complexity(i) = theoretical_complexity(xNTotal, normL, corrMod, kAdjust);
%     
% end

stream = RandStream.getDefaultStream
reset(stream, seed);

x = rand(1,5)

x = rand(1,5)
reset(stream, seed)
x = rand(1,5)

x = rand(1,5)