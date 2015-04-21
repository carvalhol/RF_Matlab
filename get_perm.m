function pVec = get_perm(pos, qmax, nStep, qmin)

        nDim = size(nStep,2);
        pVec = zeros(nDim,1);
            
        for j = 1: nDim
            seedStep = prod(nStep(j+1:end));
            if (j == nDim)
                seedStep = 1;
            end
            i = mod(floor((pos-0.9)/seedStep), nStep(j))+1;
            pVec(j) = (i-1)*(qmax(j)-qmin(j))/(nStep(j)-1)...
                       + qmin(j);
        end
end

