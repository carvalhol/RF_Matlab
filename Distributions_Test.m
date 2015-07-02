clear all
close all
clc

nIter = 100
sizeMax = 7
margifirst = 'L'; %L for Lognormal; G for Gaussian
fieldAvg = 48.6E6;
fieldVar = 10*(fieldAvg)^2;
normCut = 1E20;
normalVar = log(1 + fieldVar/(fieldAvg^2));
normalAvg = log(fieldAvg) - normalVar/2;
nBins = 100;

varGlob = zeros(sizeMax, nIter);
avgGlob = zeros(sizeMax, nIter);

for i =1:sizeMax
    
    size = 10^i;
    
    for it = 1:nIter
        
        %Std_gaussian Variable
        a = randn(size, 1);
        average = mean(a);
        stdDev = std(a);
        %hist(a,nBins)
        
        %Normalization
        % if(strcmp(margifirst, 'L'))
        %     b = sum((a>normCut))
        %     a(a>normCut) = 0;
        % end
        % a = (a - average)/stdDev;
        
        %average = mean(a)
        %variance = var(a)
        
        %Transformation
        if(strcmp(margifirst, 'L'))
            a = a * sqrt(normalVar) + normalAvg;
            a = exp(a);
        elseif(strcmp(margifirst, 'G'))
            a = a * sqrt(fieldVar) + fieldAvg;
        end
        
        average = mean(a);
        variance = var(a);        
        
        varGlob(i,it) = variance;
        avgGlob(i,it) = average;
        
    end
end

avgPercent = avgGlob/fieldAvg;
varPercent = varGlob/fieldVar;

maxErrorAvg = max(abs(avgPercent),[],2)
maxErrorVar = max(abs(varPercent),[],2)

%Method inside lognrnd
r = exp(randn(size,1) .* sqrt(normalVar) + normalAvg);
average_r = mean(r);
variance_r = var(r);
per_average_r = average_r/fieldAvg
per_variance_r = variance_r/fieldVar

%lognrnd direct utilisation
l = lognrnd(normalAvg,sqrt(normalVar),[size,1]);
average_l = mean(l);
variance_l = var(l);
per_average_l = average_l/fieldAvg
per_variance_l = variance_l/fieldVar