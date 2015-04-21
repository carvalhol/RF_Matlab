clear all
close all
clc

%% Variables
%USER
corrMod = 'G'; %G for Gaussian
nDim = 2; %Number of Dimensions
L_Glob = 2; %Starts in 0 and is a square (corrL = 1)
Nmc = 1; %Number of events
method = 'S'; %S for Shinozuka
parallelLevel = 2; % >0, will define the number of processors (fake)
independentFields = true;
seed = 1;
overLap = 0.2; %overLap in [corrL], should be an even multiple of xStep

%FIXED
corrL = 1;
xStep = corrL/10;
kAdjust = 10; %How many times the minimal k
periodMult = 1.1;
nb_procs = nDim^parallelLevel;

%OTHERS
count = 0;
xOr = zeros(nDim, 1);
stream = RandStream.getDefaultStream;

%%

figure(1)
hold on
hold all

for rank = 0:nb_procs-1
    %% Defining Local Extremes
    xProcDelta = L_Glob/parallelLevel;
    Lref = L_Glob/parallelLevel;
    L_Loc = zeros(nDim, 1);
    
    xOrVec  = 0:Lref:(L_Glob-Lref);
    
    for j=1:nDim
        seedStep = size(xOrVec,2)^(nDim-j);
        %i = floor((pos-0.9)/seedStep)+1;
        i = floor((rank+1)/seedStep)+1;
        if(mod(rank+1, seedStep) == 0)
            i = i-1;
        end
        i = mod(i-1, size(xOrVec,2))+1;
        xOr(j) = xOrVec(i);
    end      
    
    rank
    xOr
%          if (nb_procs == baseStep**nDim) then
%             write(get_fileId(),*) "Exact Division"
%             basicStep  = nint(dble(nb_procs)**(1.0d0/nDim))
%             xProcDelta = (xMaxGlob-xMinGlob)/basicStep
% 
%             !        if(rang == testRang) write(*,*) "nDim = ", nDim
%             !        if(rang == testRang) write(*,*) "basicStep = ", basicStep
%             !        if(rang == testRang) write(*,*) "xProcDelta = ", xProcDelta
%             !        if(rang == testRang) write(*,*) "nb_procs = ", nb_procs
%             !        if(rang == testRang) write(*,*) "dble(nb_procs)**(1/nDim) = ", dble(nb_procs)**(1/nDim)
%             !        if(rang == testRang) write(*,*) "nint(dble(nb_procs)**(1/nDim)) = ", nint(dble(nb_procs)**(1/nDim))
% 
%             do j = 1, nDim
%                 seedStep = basicStep**(nDim-j);
%                 i = cyclicMod(int(rang/seedStep) + 1, basicStep)
%                 !if(rang == testRang) write(*,*) "i = ", i, " seedStep = ", seedStep, "seedStep", seedStep
%                 xMin(j) = (dble(i-1))*xProcDelta(j);
%                 xMax(j) = xMin(j) + xProcDelta(j)
%                 !if(rang == testRang) write(*,*) "AFTER xMin ", j, " = ", xMin(j)
%             end do      
    %% Defining xPoints
     
    %xVec  = 0:xStep:L;
    xBottom = zeros(nDim, 1);
    xTop = zeros(nDim, 1);
    xNStep = zeros(nDim, 1);
    %xPoints = zeros(size(xVec,2)^nDim,nDim);

    if(independentFields)
        for j=1:nDim
            if(xOr(j)+Lref ~= L_Glob)
                L_Loc(j) = Lref + corrL*overLap;
            else
                L_Loc(j) = Lref;
            end
            if(xOr(j) ~= 0)
                xOr(j) = xOr(j) - corrL*overLap/2;
            end
        end
    end  
    
    for i = 1:nDim
        norm = L_Glob/xStep;
        
        xBottom(i) = xStep * ceil(xOr(i)/xStep);
        xTop(i)    = xStep * floor((xOr(i)+L_Loc(i))/xStep);
        
        
        
        %         elseif(xTop(i) == xOr(i)+L_Loc(i))
        %             if(xTop(i) ~= L_Glob)
        %                 xTop(i) = xTop(i) - xStep;
        %             end
    end
    
    xNStep = 1 + round((xTop-xBottom)./xStep);
    xNTotal = prod(xNStep);
    xPoints = zeros(xNTotal,nDim);
    
    for pos = 1:size(xPoints,1)
        xPoints(pos,:) = get_Permutation(pos, xTop + xBottom, xNStep', xBottom);
    end
    
    xOr = xBottom;
    L_Loc = xTop - xBottom;
   
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
    
    if(independentFields)
        LforK = L_Loc;
    else
        LforK = L_Glob;
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
                kPoints(pos,:) = get_Permutation(pos, kMax, kNStep', kMin);
            end
            
%             for pos = 1:size(kPoints,1)
%                 for j=1:nDim
%                     seedStep = size(kVec,2)^(nDim-j);
%                     %i = floor((pos-0.9)/seedStep)+1;
%                     i = floor(pos/seedStep)+1;
%                     if(mod(pos, seedStep) == 0)
%                         i = i-1;
%                     end
%                     i = mod(i-1, size(kVec,2))+1;
%                     kPoints(pos,j) = kVec(1,i);
%                 end
%             end
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
    
%     %% Ploting Results
%     switch nDim
%         case 1
% 
%             maxPlot = size(y,2);
%             if(maxPlot > 10)
%                 maxPlot = 10;
%             end
%             for event=1:maxPlot
%                 plot(xPoints, y(:,event))
%             end
%             %hold off
%         case 2
%             %figure(1)
%             %hold on
%             maxPlot = 1;
%             for event=1:maxPlot
%                 plot2D = plot3(xPoints(:,1), xPoints(:,2), y(:,event),'o');
%                 tri2D = delaunay(xPoints(:,1), xPoints(:,2));
%                 trisurf(tri2D,xPoints(:,1), xPoints(:,2), y(:,event))
%                 delete(plot2D)
%             end
%             %hold off
%         case 3
%             %figure(1)
%             %hold on
%             maxPlot = 1;
%             for event=1:maxPlot
%                 %plot3D = plot3(xPoints(:,1), xPoints(:,2), xPoints(:,3),'o');
%                 %tri3D  = delaunay(xPoints(:,1), xPoints(:,2), xPoints(:,3));
%                 %trisurf(tri3D,xPoints(:,1), xPoints(:,2), xPoints(:,3), y(:,event))
%                 scatter3(xPoints(:,1), xPoints(:,2), xPoints(:,3),15, y(:,event),'filled')
%             end
%             %hold off
%             
%         otherwise
%             error('Dimension not accepted in plot')
%     end
    
end

hold off

% %% Statistics
% 
% average = mean(mean(y))
% stdDev = mean(std(y))
% 1/sqrt(size(y,1))

