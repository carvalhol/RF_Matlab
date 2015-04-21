clear all
close all
clc

%% Variables
%USER
corrMod = 'G'; %G for Gaussian
nDim = 2; %Number of Dimensions
L = 100; %Starts in 0 and is a square (corrL = 1)
Nmc = 1; %Number of events
method = 'S'; %S for Shinozuka

%FIXED
corrL = 1;
xStep = corrL/10;
kAdjust = 10; %How many times the minimal k
periodMult = 1.1;

%OTHERS
xVec  = 0:xStep:L;
xPoints = zeros(size(xVec,2)^nDim,nDim);
count = 0;

%% Defining xPoints            
for pos = 1:size(xPoints,1)
    for j=1:nDim
        seedStep = size(xVec,2)^(nDim-j);
        %i = floor((pos-0.9)/seedStep)+1;
        i = floor(pos/seedStep)+1;
        if(mod(pos, seedStep) == 0)
            i = i-1;
        end
        i = mod(i-1, size(xVec,2))+1;
        xPoints(pos,j) = xVec(1,i);
    end
end

%% Defining kMax

integralExigence = 0.01;
kmax = 1;
proportion = 1;

switch corrMod
    case 'G'
        ref = pi^nDim;
        
        switch nDim                                   
            case 1
                kMax = 6.457; %integralExigence = 0.01
%                 fun = @(x) exp(-x.^2/(4*pi));                             
%                 while(proportion > integralExigence)
%                     kMax = kmax+0.001;
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

% est_kPointsSize = (1 + kMax/(kMax/ceil(kMax/(2*pi/L))))^nDim;
% est_kPointsSize = (1 + (ceil(kMax/(2*pi/L))))^nDim;
% est_kPointsSize = (1 + (ceil(L*kMax/(2*pi))))^nDim;
est_kPointsSize = (1 + kAdjust*(ceil(L*kMax/(2*pi))))^nDim;

%% Defining kPoints
switch method
    case 'S'
        kStepMax = 2*pi/(periodMult*L);
        kStep = kMax/ceil(kMax/kStepMax);
        kVec  = 0:kStep/kAdjust:kMax;
        kPoints = zeros(size(kVec,2)^nDim,nDim);
        for pos = 1:size(kPoints,1)
            for j=1:nDim
                seedStep = size(kVec,2)^(nDim-j);
                %i = floor((pos-0.9)/seedStep)+1;
                i = floor(pos/seedStep)+1;
                if(mod(pos, seedStep) == 0)
                    i = i-1;
                end
                i = mod(i-1, size(kVec,2))+1;
                kPoints(pos,j) = kVec(1,i);
            end
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
switch method
    case 'S'
        y = zeros(size(xPoints,1),Nmc);
        ampMult = 2*sqrt((kStep/(2*pi))^(nDim));
        %         2.0d0*sqrt(product(kDelta)/((2.0d0*PI)**(dble(RDF%nDim))))
        for event = 1:Nmc            
            %phiK = zeros(1, size(kPoints,1));
            %phiK = unifrnd(0,2*pi,1, size(kPoints,1));
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
        figure(1)
        hold on
        hold all
        maxPlot = size(y,2);
        if(maxPlot > 10)
            maxPlot = 10;
        end
        for event=1:maxPlot
            plot(xPoints, y(:,event))
        end
        hold off
    case 2
        figure(1)
        hold on
        maxPlot = 1;
        for event=1:maxPlot
            plot2D = plot3(xPoints(:,1), xPoints(:,2), y(:,event),'o');
            tri2D = delaunay(xPoints(:,1), xPoints(:,2));
            trisurf(tri2D,xPoints(:,1), xPoints(:,2), y(:,event))
            delete(plot2D)
        end
        hold off
    case 3
        figure(1)
        hold on
        maxPlot = 1;
        for event=1:maxPlot
            %plot3D = plot3(xPoints(:,1), xPoints(:,2), xPoints(:,3),'o');
            %tri3D  = delaunay(xPoints(:,1), xPoints(:,2), xPoints(:,3));
            %trisurf(tri3D,xPoints(:,1), xPoints(:,2), xPoints(:,3), y(:,event))
            scatter3(xPoints(:,1), xPoints(:,2), xPoints(:,3),15, y(:,event),'filled')
        end
        hold off
        
    otherwise
        error('Dimension not accepted in plot')
end

%% Statistics

average = mean(mean(y))
stdDev = mean(std(y))
1/sqrt(size(y,1))

