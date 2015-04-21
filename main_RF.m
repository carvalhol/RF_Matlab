clear all
close all
clc

%% Variables
%USER
corrMod = 'G'; %G for Gaussian
nDim = 2; %Number of Dimensions
L_Glob = 4; %Starts in 0 and is a square (corrL = 1)
Nmc = 1; %Number of events
method = 'S'; %S for Shinozuka
parallelLevel = 2; % >0, will define the number of processors (fake)
independentFields = true;
seed = 1;
overLap = 0.4; %overLap in [corrL], should be an even multiple of xStep

%FIXED
corrL = 1;
xStep = corrL/10;
kAdjust = 10; %How many times the minimal k
periodMult = 1.1;
nb_procs = parallelLevel^nDim;

%OTHERS
count = 0;
xMin = zeros(nDim, 1);
xMax = zeros(nDim, 1);
stream = RandStream.getDefaultStream;

%%

figure(1)
hold on
hold all

for rank = 0:nb_procs-1
    %% Defining Local Extremes
    xProcDelta = L_Glob/parallelLevel;
    L_Loc = zeros(nDim, 1);
    
    xMinSeed  = 0:xProcDelta:(L_Glob-xProcDelta);
    
    for j=1:nDim
        seedStep = size(xMinSeed,2)^(nDim-j);
        %i = floor((pos-0.9)/seedStep)+1;
        i = floor((rank+1)/seedStep)+1;
        if(mod(rank+1, seedStep) == 0)
            i = i-1;
        end
        i = mod(i-1, size(xMinSeed,2))+1;
        xMin(j) = xMinSeed(i);
        xMax(j) = xMin(j) + xProcDelta;
    end      
    
%     rank
%     xMin
%     xMax
 
    %% Defining xPoints borders
    
    xMaxRef = xMax(:);

    if(independentFields)
        for j=1:nDim
            if(xMax(j) ~= L_Glob)
               xMax(j) = xMax(j) + corrL*overLap/2;
            end
            if(xMin(j) ~= 0)
                xMin(j) = xMin(j) - corrL*overLap/2;
            end
        end
    end  
    
    %Rounding to nearest point on the grid
    for j = 1:nDim
        
        xMin(j) = xStep * ceil(xMin(j)/xStep);
        xMax(j) = xStep * floor((xMax(j))/xStep);
                       
        if(xMax(j) == xMaxRef(j))
            if(xMaxRef(j) ~= L_Glob)
                xMax(j) = xMax(j) - xStep;
            end
        end
    end
    
    L_Loc = xMax(:) - xMin(:);
    
%     rank
%     xMin
%     xMax
%     L_Loc
    
    %% Defining xPoints
    xNStep = zeros(nDim, 1);
    xNStep = 1 + round((xMax-xMin)./xStep);
    xNTotal = prod(xNStep);
    xPoints = zeros(xNTotal,nDim);
    
    for pos = 1:size(xPoints,1)
        xPoints(pos,:) = get_perm(pos, xMax, xNStep', xMin);
    end
    
    
   
    %% Defining kMax
    
    integralExigence = 0.01;
    kMax = zeros(nDim, 1);
    kMin = zeros(nDim, 1);
    proportion = 1;
    
    switch corrMod
        case 'G'
            ref = pi^nDim;
            
            switch nDim
                case 1
                    kMax(:) = 6.457; %integralExigence = 0.01
                    %                 fun = @(x) exp(-x.^2/(4*pi));
                    %                 while(proportion > integralExigence)
                    %                     kMax = kmax+0.001;
                    %                     trial = quad(fun,0,kMax);
                    %                     proportion = 1 - trial/ref;
                    %                 end
                case 2
                    kMax(:) = 7.035; %integralExigence = 0.01
                    %                 fun = @(x,y) exp(-(x.^2+y.^2)/(4*pi));
                    %                 while(proportion > integralExigence)
                    %                     kMax = kMax+0.001;
                    %                     trial = quad2d(fun,0,kMax,0,kMax);
                    %                     proportion = 1 - trial/ref;
                    %                 end
                case 3
                    kMax(:) = 7.355; %integralExigence = 0.01
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
    
    % est_kPointsSize = (1 + kMax/(kMax/ceil(kMax/(2*pi/L))))^nDim;
    % est_kPointsSize = (1 + (ceil(kMax/(2*pi/L))))^nDim;
    % est_kPointsSize = (1 + (ceil(L*kMax/(2*pi))))^nDim;
    
    
    %% Defining kPoints
    
    LforK = L_Loc;
    if(~independentFields)
        LforK(:) = L_Glob;
    end
    
    est_kPointsSize = prod(1 + kAdjust*(ceil(LforK.*kMax/(2*pi))));%^nDim;
    
    switch method
        case 'S'
            kStepMax = 2*pi./(periodMult*LforK);
            kStep = kMax./ceil(kMax./kStepMax);
            kNStep = kMax./kStep;
            %kVec  = 0:kStep/kAdjust:kMax;
            kPoints = zeros(prod(kNStep),nDim);
            
            for pos = 1:size(kPoints,1)
                kPoints(pos,:) = get_perm(pos, kMax, kNStep', kMin);
            end
            
        otherwise
            error('method not implemented')
    end
    
    %% Defining Spectrum
    switch corrMod
        case 'G'
            Sk = exp(-sum(kPoints.^2, 2)/(4*pi));
        otherwise
            error('corMod not implemented')
    end
    
    %% Calculating Random Field
    if(~independentFields)
        reset(stream, seed);
    end    
    
    switch method
        case 'S'
            y = zeros(size(xPoints,1),Nmc);
            ampMult = prod(2*sqrt((kStep./(2*pi))));%^(nDim));
            for event = 1:Nmc
                phiK = 2*pi*rand(1, size(kPoints,1));
                y(:,event) = cos(xPoints*kPoints' + repmat(phiK,size(xPoints,1),1))*sqrt(Sk);
            end
            
            y = ampMult*y;
            
        otherwise
            error('method not implemented')
    end
    
    %% Ploting Results
    switch nDim
        case 1

            maxPlot = size(y,2);
            if(maxPlot > 10)
                maxPlot = 10;
            end
            for event=1:maxPlot
                plot(xPoints, y(:,event))
            end
            %hold off
        case 2
            %figure(1)
            %hold on
            maxPlot = 1;
            for event=1:maxPlot
                plot2D = plot3(xPoints(:,1), xPoints(:,2), y(:,event),'o');
                tri2D = delaunay(xPoints(:,1), xPoints(:,2));
                trisurf(tri2D,xPoints(:,1), xPoints(:,2), y(:,event))
                delete(plot2D)
            end
            %hold off
        case 3
            %figure(1)
            %hold on
            maxPlot = 1;
            for event=1:maxPlot
                %plot3D = plot3(xPoints(:,1), xPoints(:,2), xPoints(:,3),'o');
                %tri3D  = delaunay(xPoints(:,1), xPoints(:,2), xPoints(:,3));
                %trisurf(tri3D,xPoints(:,1), xPoints(:,2), xPoints(:,3), y(:,event))
                scatter3(xPoints(:,1), xPoints(:,2), xPoints(:,3),15, y(:,event),'filled')
            end
            %hold off
            
        otherwise
            error('Dimension not accepted in plot')
    end
    
end

hold off

% %% Statistics
% 
% average = mean(mean(y))
% stdDev = mean(std(y))
% 1/sqrt(size(y,1))

