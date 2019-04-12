%% SOME CONSTANTS
clear;
clc; %clear any previous workspaces

windowLength = 0.03; %length of window is seconds (30ms)
overlap = 0.5; %percentage of window overlap
Po = 1.2; %kg/m^3
c = 343; %m/s
Zo = Po*c; % = 411.6 Pa.s/m

% SET THESE TO 'y' IF WANT TO CHECK, SET AS 'n' TO NOT TEST AND SAVE TIME
qualityCheck = 'y';
localisationCheck = 'y';
diffusionCheck = 'y';

%% READING THROUGH AMBISONIC FOLDERS
% BASED OFF "https://uk.mathworks.com/matlabcentral/answers/uploaded_files/30598/recurse_subfolders.m"

% Define a starting folder.
topLevelFolder = 'C:\Users\david\Documents\Result';

allSubFolders = genpath(topLevelFolder);
% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {}; % create folder name array
while true
	[singleSubFolder, remain] = strtok(remain, ';');
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames);  

%% PROCESSING 

for Folder=2:numberOfFolders % first folder in list is top level folder, ignore
    cd (listOfFolderNames{Folder}) % changes folder for each loop
    
    mylist = dir('*.wav'); %make list of wav files in this folder
    listlength = length(mylist); %find the length of this list 
    mylist = struct2cell(mylist); % need cell, not struct to read eachfile with audioread 
    
    fileAziAngle = zeros(listlength,1);
    fileSSIM = zeros(listlength,1);
    fileDiffusion = zeros(listlength,1);
    
    %% Read In File
    for listIndex= 1:listlength %for the amount of files in this folder
        % Ambisonic file read in
        fn(1,1)= mylist(1,listIndex); 
        [x1,fs]=audioread(fn{1,1}); %read current file, store vals in x, fs = sampling frequency
        sigLength = length(x1); %store the length of this wav file 
        % which order Ambisonics is it?
        width = size(x1,2); % finds # channels and therefore order
        
        origambi = mylist(1,1);
        [orig, fs] = audioread(origambi{1,1});
        
        if width == 10
            width = 9;
        end
        
        %% SET-UP CALCS FOR FILE
        % general
        runtime= round((length(x1)-1)/fs); % find overal time of clip in secs
        numsamp = floor(windowLength*fs); %calculate #samples per window
        update = floor(fs*overlap); %find #samples for this overlap time 

        arraySize = floor(length(x1)/update); %how many windows for this clip?
        maxXval = arraySize*update; % will truncate some of the signal to round value
        xpoints = (numsamp/2:update:maxXval)/fs; % times of center of each window for plotting
        
        %localisation
        aziCoeffAvg = zeros(arraySize,width);
        aziFuMaCoeffs = [1;1;2/sqrt(3);2/sqrt(3);sqrt(45/32); sqrt(45/32);sqrt(8/5);sqrt(8/5)];
        aziCols = [2,3,8,9,11,12,15,16]; %which channels/columns to use for azimuth

        aziWindowAngle = zeros(arraySize,1); %stores average angle per frame for all frames

        omniFuMaCoeff = [1/sqrt(2)];
        elevFuMaCoeffs = [1;1;2/sqrt(3);2/sqrt(3);1; 3/sqrt(5); 3/sqrt(5);];

        % timbral quality
        SSIMcoeffs = zeros(arraySize,1);
        
        % Diffusion arrays
        %IACC = zeros(arraySize,1);
        diffuseness = zeros(arraySize,1);
        
        % for processing
        pos =1; %index value, update value
        window = 1; % index for windows
        
        switch width %depending on order do this
        %% FIRST ORDER FILE        
            case 4 % for a first order file
                aziWidth = 2;
                x = [x1(:,1),x1(:,4),x1(:,2),x1(:,3)]; %change from ambiX order to FuMa order
                while (pos+numsamp) <= length(x) % to the end of the sampled data
                    y=x(pos:pos+numsamp-1, :); % y has the x sampled data of this window
                    t =(0:length(y)-1)/fs; % finds the time each sample is at 
                    origy = orig(pos:pos+numsamp-1, :);
                    
                    %% LOCALISATION ACCURACY
                    if localisationCheck == 'y'
                        for check=1:length(y) %omni value may be negative, flip for calcs
                           if y(check,1) < 0 
                               y(check, 2) = y(check,2)*(-1);
                               y(check, 3) = y(check,3)*(-1);
                               y(check, 4) = y(check,4)*(-1);
                           end
                        end
                
                        for j=1:length(y)
                            for i=1:aziWidth
                                aziCoeffAvg(j,i) = y(j,aziCols(1,i))*aziFuMaCoeffs(i,1);
                                %multiply each channel value by the respective FuMa coefficient
                            end 
                        end
                        
                        aziCoeffAvg(1,1) = sum(aziCoeffAvg(:,1))/length(aziCoeffAvg);
                        aziCoeffAvg(1,2) = sum(aziCoeffAvg(:,2))/length(aziCoeffAvg);
                        
                        i=1;
                        for t=0:0.01:2*pi
                            vals(i,1) = (aziCoeffsAvg(1,1)*cos(t)) + (aziCoeffsAvg(1,2)*sin(t));
                            xvals(i,1) = t;
                            i = i+1;
                        end

                        maxrad = max(vals);
                        for i=1:length(vals)
                            if maxrad == vals(i,1)
                               avgAngle = 180*(xvals(i,1))/pi;
                            end
                        end
                        
                        aziWindowAngle(window,1) = avgAngle;
                    end
                      
                    %% TIMBRAL QUALITY
                    
                    if qualityCheck == 'y'
                        WchannelOrig = origy(:,1);
                        WchannelTest = y(:,1);
                        spgrambw(WchannelOrig,fs,'Jwta');
                        saveas(gcf,'ref.png');
                        spgrambw(WchannelTest,fs,'Jwta');
                        saveas(gcf,'test.png');
                        close all;
                        
                        compref = imread('ref.png');
                        comptest = imread('test.png');
                        ssimval = ssim(compref,comptest);
                        delete *.png;
                        
                        SSIMcoeffs(window,1) = ssimval;
                    end
                    
                    %% DIFFUSION RATING
                    
                    if diffusionCheck == 'y'
                        % I = pu, where p = W channel, U = velocity vector
                        p = y(:,1); %take W channel only
                        U = zeros(length(y),1); 

                        %velocity vectors of each channel
                        uy = (1/(Zo*sqrt(2)))*(y(:,2));
                        uz = (1/(Zo*sqrt(2)))*(y(:,3));
                        ux = (1/(Zo*sqrt(2)))*(y(:,4));

                        %find combined velocity vector
                        for i =1:length(y)
                            U(i,1) = ux(i,1)+ uy(i,1)+ uz(i,1);
                        end

                        I = p.*U;
                        
                        % energy vector
                        E1 = (y(:,1).^2);
                        E2 = (y(:,2).^2);
                        E3 = (y(:,3).^2);
                        E4 = (y(:,4).^2);
                        expE1 = sum(E1)/length(E1);
                        expE2 = sum(E2)/length(E2);
                        expE3 = sum(E3)/length(E3);
                        expE4 = sum(E4)/length(E4);
                        E = Po*(1/2)*(1/(Zo)^2)*(expE1 + (expE2+expE3+expE4)/2);
                        %E = Po*(1/2)*((expE1/(Zo^2))*avgU);
                        
                        ExpI = sum(I)/length(I);

                        diffuseness(window, 1) = 1 - (abs(ExpI)/(c*E));
                    end

                    % update
                    pos=pos+update; % This gives an overlap window of Y%
                    window = window +1; %update window counter
                end %end of window
               
                
        %% SECOND ORDER FILE
            case 9 % for a second order file
                aziWidth = 4;
                x = [x1(:,1),x1(:,4),x1(:,2),x1(:,3),x1(:,7),x1(:,8),x1(:,6),x1(:,9),x1(:,4)]; %change from ambiX order to FuMa order
                while (pos+numsamp) <= length(x) % to the end of the sampled data
                    y=x(pos:pos+numsamp-1, :); % y has the x sampled data of this window
                    t =(0:length(y)-1)/fs; % finds the time each sample is at 
                    origy = orig(pos:pos+numsamp-1, :);
                    
                    %% LOCALISATION ACCURACY
                    if localisationCheck == 'y'
                        for check=1:length(y) %omni value may be negative, flip for calcs
                           if y(check,1) < 0 
                               y(check, 2) = y(check,2)*(-1);
                               y(check, 3) = y(check,3)*(-1);
                               y(check, 4) = y(check,4)*(-1);
                               y(check, 5) = y(check,5)*(-1);
                               y(check, 6) = y(check,6)*(-1);
                               y(check, 7) = y(check,7)*(-1);
                               y(check, 8) = y(check,8)*(-1);
                               y(check, 9) = y(check,9)*(-1);
                           end
                        end

                        for j=1:length(y)
                            for i=1:aziWidth
                                aziCoeffAvg(j,i) = y(j,aziCols(1,i))*aziFuMaCoeffs(i,1);
                                %multiply each channel value by the respective FuMa coefficient
                            end 
                        end
                        
                        aziCoeffAvg(1,1) = sum(aziCoeffAvg(:,1))/length(aziCoeffAvg);
                        aziCoeffAvg(1,2) = sum(aziCoeffAvg(:,2))/length(aziCoeffAvg);
                        aziCoeffAvg(1,3) = sum(aziCoeffAvg(:,3))/length(aziCoeffAvg);
                        aziCoeffAvg(1,4) = sum(aziCoeffAvg(:,4))/length(aziCoeffAvg);
                        
                        i = 1;   
                        for t=0:0.01:2*pi
                            vals(i,1) = (aziCoeffAvg(1,1)*cos(t)) + (aziCoeffAvg(1,2)*sin(t)) + (aziCoeffAvg(1,3)*(sqrt(3)/2)*(cos(2*t))) + (aziCoeffAvg(1,4)*(sqrt(3)/2)*(sin(2*t)));
                            xvals(i,1) = t;
                            i = i+1;
                        end

                        maxrad = max(vals);
                        for i=1:length(vals)
                            if maxrad == vals(i,1)
                               aziAngle = 180*(xvals(i,1))/pi;
                            end
                        end
                        
                        aziWindowAngle(window,1) = aziAngle;
                    end
                       
                    %% TIMBRAL QUALITY
                    
                    if qualityCheck == 'y'
                        WchannelOrig = origy(:,1);
                        WchannelTest = y(:,1);
                        spgrambw(WchannelOrig,fs,'Jwta');
                        saveas(gcf,'ref.png');
                        spgrambw(WchannelTest,fs,'Jwta');
                        saveas(gcf,'test.png');
                        close all;
                        
                        compref = imread('ref.png');
                        comptest = imread('test.png');
                        ssimval = ssim(compref,comptest);
                        delete *.png;
                        
                        SSIMcoeffs(window,1) = ssimval;
                    end
                    
                    %% DIFFUSION RATING
                    
                    if diffusionCheck == 'y'
                        % I = pu, where p = W channel, U = velocity vector
                        p = y(:,1); %take W channel only
                        U = zeros(length(y),1); 

                        %velocity vectors of each channel
                        uy = (1/(Zo*sqrt(2)))*(y(:,2));
                        uz = (1/(Zo*sqrt(2)))*(y(:,3));
                        ux = (1/(Zo*sqrt(2)))*(y(:,4));
                        uu = (1/(Zo*sqrt(2)))*(y(:,5));
                        uv = (1/(Zo*sqrt(2)))*(y(:,6));
                        us = (1/(Zo*sqrt(2)))*(y(:,7));
                        ut = (1/(Zo*sqrt(2)))*(y(:,8));
                        ur = (1/(Zo*sqrt(2)))*(y(:,9));

                        %find combined velocity vector
                        for i =1:length(y)
                            U(i,1) = ux(i,1)+ uy(i,1)+ uz(i,1) + uu(i,1) + uv(i,1) + us(i,1) + ut(i,1) + ur(i,1);
                        end

                        I = p.*U;
                        
                        E1 = (y(:,1).^2);
                        E2 = (y(:,2).^2);
                        E3 = (y(:,3).^2);
                        E4 = (y(:,4).^2);
                        E5 = (y(:,5).^2);
                        E6 = (y(:,6).^2);
                        E7 = (y(:,7).^2);
                        E8 = (y(:,8).^2);
                        E9 = (y(:,9).^2);
                        expE1 = sum(E1)/length(E1);
                        expE2 = sum(E2)/length(E2);
                        expE3 = sum(E3)/length(E3);
                        expE4 = sum(E4)/length(E4);
                        expE5 = sum(E5)/length(E5);
                        expE6 = sum(E6)/length(E6);
                        expE7 = sum(E7)/length(E7);
                        expE8 = sum(E8)/length(E8);
                        expE9 = sum(E9)/length(E9);
                        E = Po*(1/2)*(1/(Zo^2))*(expE1 + (expE2+expE3+expE4+expE5+expE6+expE7+expE8+expE9)/2);

                        ExpI = sum(I)/length(I);

                        diffuseness(window, 1) = 1 - (abs(ExpI)/(c*E));
                    end
                    
                    % update
                    pos=pos+update; % This gives an overlap window of Y%
                    window = window +1; %update window counter 
                end % end of window
                
                
                
        %% THIRD ORDER FILE
            case 16 % for a third order file
                aziWidth = 8;
                x = [x1(:,1),x1(:,4),x1(:,2),x1(:,3),x1(:,7),x1(:,8),x1(:,6),x1(:,9),x1(:,4),x1(:,13),x1(:,14),x1(:,12),x1(:,15),x1(:,11),x1(:,16),x1(:,10)]; %change from ambiX order to FuMa order
                while (pos+numsamp) <= length(x) % to the end of the sampled data
                    y=x(pos:pos+numsamp-1, :); % y has the x sampled data of this window
                    t =(0:length(y)-1)/fs; % finds the time each sample is at 
                    origy = orig(pos:pos+numsamp-1, :);
                    
                    for check=1:length(y) %omni value may be negative, flip for calcs
                       if y(check,1) < 0 
                           y(check, 2) = y(check,2)*(-1);
                           y(check, 3) = y(check,3)*(-1);
                           y(check, 4) = y(check,4)*(-1);
                           y(check, 5) = y(check,5)*(-1);
                           y(check, 6) = y(check,6)*(-1);
                           y(check, 7) = y(check,7)*(-1);
                           y(check, 8) = y(check,8)*(-1);
                           y(check, 9) = y(check,9)*(-1);
                           y(check, 10) = y(check,10)*(-1);
                           y(check, 11) = y(check,11)*(-1);
                           y(check, 12) = y(check,12)*(-1);
                           y(check, 13) = y(check,13)*(-1);
                           y(check, 14) = y(check,14)*(-1);
                           y(check, 15) = y(check,15)*(-1);
                           y(check, 16) = y(check,16)*(-1);
                       end
                    end
                
                    %% LOCALISATION ACCURACY
                    if localisationCheck == 'y'
                        for j=1:length(y)
                            for i=1:aziWidth
                                aziCoeffAvg(j,i) = y(j,aziCols(1,i))*aziFuMaCoeffs(i,1);
                                %multiply each channel value by the respective FuMa coefficient
                            end 
                        end
                        
                        aziCoeffAvg(1,1) = sum(aziCoeffAvg(:,1))/length(aziCoeffAvg);
                        aziCoeffAvg(1,2) = sum(aziCoeffAvg(:,2))/length(aziCoeffAvg);
                        aziCoeffAvg(1,3) = sum(aziCoeffAvg(:,3))/length(aziCoeffAvg);
                        aziCoeffAvg(1,4) = sum(aziCoeffAvg(:,4))/length(aziCoeffAvg);
                        aziCoeffAvg(1,5) = sum(aziCoeffAvg(:,5))/length(aziCoeffAvg);
                        aziCoeffAvg(1,6) = sum(aziCoeffAvg(:,6))/length(aziCoeffAvg);
                        aziCoeffAvg(1,7) = sum(aziCoeffAvg(:,7))/length(aziCoeffAvg);
                        aziCoeffAvg(1,8) = sum(aziCoeffAvg(:,8))/length(aziCoeffAvg);
                        
                        i=1;
                        for t=0:0.01:2*pi
                            vals(i,1) = (aziCoeffAvg(1,1)*cos(t)) + (aziCoeffAvg(1,2)*sin(t)) + (aziCoeffAvg(1,3)*(sqrt(3)/2)*(cos(2*t))) + (aziCoeffAvg(1,4)*(sqrt(3)/2)*(sin(2*t))) + (aziCoeffAvg(1,5)*(-sqrt(3)/8)*(cos(t))) + (aziCoeffAvg(1,6)*(-sqrt(3)/8)*(sin(t))) + (aziCoeffAvg(1,7)*(sqrt(5)/8)*(cos(3*t))) + (aziCoeffAvg(1,8)*(sqrt(5)/8)*(sin(3*t)));
                            xvals(i,1) = t;
                            i = i+1;
                        end

                        maxrad = max(vals);
                        for i=1:length(vals)
                            if maxrad == vals(i,1)
                               aziAngle = 180*(xvals(i,1))/pi;
                            end
                        end

                        aziWindowAngle(window,1) = aziAngle;
                    end
                    
                    %% TIMBRAL QUALITY
                    
                    if qualityCheck == 'y'
                        WchannelOrig = origy(:,1);
                        WchannelTest = y(:,1);
                        spgrambw(WchannelOrig,fs,'Jwta');
                        saveas(gcf,'ref.png');
                        spgrambw(WchannelTest,fs,'Jwta');
                        saveas(gcf,'test.png');
                        close all;
                        
                        compref = imread('ref.png');
                        comptest = imread('test.png');
                        ssimval = ssim(compref,comptest);
                        delete *.png;
                        
                        SSIMcoeffs(window,1) = ssimval;
                    end
                    
                    %% DIFFUSION RATING
                    if diffusionCheck == 'y'
                        % I = pu, where p = W channel, U = velocity vector
                        p = y(:,1); %take W channel only
                        U = zeros(length(y),1); 

                        %velocity vectors of each channel
                        uy = (1/(Zo*sqrt(2)))*(y(:,2));
                        uz = (1/(Zo*sqrt(2)))*(y(:,3));
                        ux = (1/(Zo*sqrt(2)))*(y(:,4));
                        uu = (1/(Zo*sqrt(2)))*(y(:,5));
                        uv = (1/(Zo*sqrt(2)))*(y(:,6));
                        us = (1/(Zo*sqrt(2)))*(y(:,7));
                        ut = (1/(Zo*sqrt(2)))*(y(:,8));
                        ur = (1/(Zo*sqrt(2)))*(y(:,9));
                        uk = (1/(Zo*sqrt(2)))*(y(:,10));
                        ul = (1/(Zo*sqrt(2)))*(y(:,11));
                        um = (1/(Zo*sqrt(2)))*(y(:,12));
                        un = (1/(Zo*sqrt(2)))*(y(:,13));
                        uo = (1/(Zo*sqrt(2)))*(y(:,14));
                        up = (1/(Zo*sqrt(2)))*(y(:,15));
                        uq = (1/(Zo*sqrt(2)))*(y(:,16));

                        %find combined velocity vector
                        for i =1:length(y)
                            U(i,1) = ux(i,1)+ uy(i,1)+ uz(i,1) + uu(i,1) + uv(i,1) + us(i,1) + ut(i,1) + ur(i,1)+ uk(i,1)+ ul(i,1) + um(i,1) + un(i,1) + uo(i,1) + up(i,1) + uq(i,1);
                        end

                        I = p.*U;

                        % energy vector
                        E1 = (y(:,1).^2);
                        E2 = (y(:,2).^2);
                        E3 = (y(:,3).^2);
                        E4 = (y(:,4).^2);
                        E5 = (y(:,5).^2);
                        E6 = (y(:,6).^2);
                        E7 = (y(:,7).^2);
                        E8 = (y(:,8).^2);
                        E9 = (y(:,9).^2);
                        E10 = (y(:,10).^2);
                        E11 = (y(:,11).^2);
                        E12 = (y(:,12).^2);
                        E13 = (y(:,13).^2);
                        E14 = (y(:,14).^2);
                        E15 = (y(:,15).^2);
                        E16 = (y(:,16).^2);
                        expE1 = sum(E1)/length(E1);
                        expE2 = sum(E2)/length(E2);
                        expE3 = sum(E3)/length(E3);
                        expE4 = sum(E4)/length(E4);
                        expE5 = sum(E5)/length(E5);
                        expE6 = sum(E6)/length(E6);
                        expE7 = sum(E7)/length(E7);
                        expE8 = sum(E8)/length(E8);
                        expE9 = sum(E9)/length(E9);
                        expE10 = sum(E10)/length(E10);
                        expE11 = sum(E11)/length(E11);
                        expE12 = sum(E12)/length(E12);
                        expE13 = sum(E13)/length(E13);
                        expE14 = sum(E14)/length(E14);
                        expE15 = sum(E15)/length(E15);
                        expE16 = sum(E16)/length(E16);
                        E = Po*(1/2)*(1/(Zo^2))*(expE1 + (expE2+expE3+expE4+expE5+expE6+expE7+expE8+expE9+expE10+expE11+expE12+expE13+expE14+expE15+expE16)/2);
                        ExpI = sum(I)/length(I);

                        diffuseness(window, 1) = 1 - (abs(ExpI)/(c*E));
                    end
                    
                    % update
                    pos=pos+update; % This gives an overlap window of Y%
                    window = window +1; %update window counter
                end % end of window
        end % end of switch 
        
        fileAziAngle(listIndex,1) = sum(aziWindowAngle)/length(aziWindowAngle);
        fileSSIM(listIndex, 1) = sum(SSIMcoeffs)/length(SSIMcoeffs);
        fileDiffusion(listIndex,1) = sum(diffuseness)/length(diffuseness);
        
    end % end of file processing
    
    allAziAngles(:,Folder-1) = fileAziAngle(:,1);
    allSSIM(:,Folder-1) = fileSSIM(:,1);
    allDiffusion(:, Folder-1) = fileDiffusion(:,1);
    
    % send end results for that original file + compressed files into overall results array
    % matrix reads: every col = different file
    % row1: original, row2: bitrate1, row3: bitrate2, row4: bitrate3 etc.
    
end % end of folder processing
