%FIGURES FIGURES FIGURES


%% GRAB FOLDERS AND FILES
%Running this first module gets the file list
clc; close all;

% Define a starting folder.
start_path = fullfile(matlabroot, '\toolbox');
if ~exist(start_path, 'dir')
	start_path = matlabroot;
end
% Ask user to confirm the folder, or change it.
uiwait(msgbox('Pick a starting folder on the next window that will come up.'));
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
	return;
end
fprintf('The top level folder is "%s".\n', topLevelFolder);

%grab parameter files
filePattern = sprintf('%s/**/*_particle_params.csv', topLevelFolder);
allFileInfo = dir(filePattern);
%add .mat files
% filePattern = sprintf('%s/**/*all_gen_params.mat', topLevelFolder);
% allFileInfo = [allFileInfo;dir(filePattern)];

% Get a cell array of strings.  
listOfFolderNames = unique({allFileInfo.folder});
numberOfFolders = length(listOfFolderNames);
fprintf('The total number of folders to look in is %d.\n', numberOfFolders);

% Get a cell array of base filename strings.  
listOfFileNames = {allFileInfo.name};
totalNumberOfFiles = length(listOfFileNames);
fprintf('The total number of files in those %d folders is %d.\n', numberOfFolders, totalNumberOfFiles);

% Process all files in those folders.
totalNumberOfFiles = length(allFileInfo);
% Now we have a list of all files, matching the pattern, in the top level folder and its subfolders.
%% Variables
R = zeros(totalNumberOfFiles*3,15); %STORE region data, 
time = zeros(totalNumberOfFiles,1); %time per video 
H = zeros(totalNumberOfFiles,1);
L = string(zeros(totalNumberOfFiles*3,1));
%% COLOR

color_path = '/Users/marysalcedo/Documents/MATLAB/Wings_Matlab_Genscode/';
addpath('/Users/marysalcedo/Documents/MATLAB/Wings_Matlab_Genscode/cbrewer');
colorB = cbrewer('div', 'RdYlBu',80); %colormap length of A
%% X Y COORDINATES and LEFT/RIGHT/ABD/THORAX/PRO LABELS
% Running this module went smoothly

j = 0; %counters for number of particles
h = 1;
count = 1; %for L/R/ABD/THX/PRO labels
for k = 1:totalNumberOfFiles
   
    % Go through all those files, grab x/y parameter file
    thisFolder = allFileInfo(k).folder;
    thisBaseFileName = allFileInfo(k).name;
    fullFileName = fullfile(thisFolder, thisBaseFileName);
    fprintf('     Processing file %d of %d : "%s".\n', k, totalNumberOfFiles, fullFileName);
    
    fileParams = dir(fullfile(thisFolder, '*xy_coord_params.mat')); %computed
    genParams = dir(fullfile(thisFolder, '*all_gen_params.mat'));
    
    xy = load(fullfile(thisFolder, fileParams.name));
    xy_coord = xy.xy_coord;
    
    gen = load(fullfile(thisFolder, genParams.name)); %get wing region values
    gen = gen.A;
    a = find(~cellfun(@isempty,gen)); %find when parameters are not empty
    reg = gen{1,a}{2,1}; %gives region as a char value
    reg = str2double(reg);
    
    
    xy_coord = cellstr(xy_coord{1,1}); %convert to cell array
    [m,~] = size(xy_coord);
    j = j+m; %update j

    for g = 1:m
        a = xy_coord(g); %take 1st val in cell array
        
        c = find(a{1,1} == ','); %use char values
        
        x(g) = str2double(a{1,1}(1:c(1)-1)); %find first comma, x values
        y(g) = str2double(a{1,1}(c(1)+1:end));   %y values
        
        % Left wings, Right wings, ABD, Heart, Pronotum, Thorax
        %X and Y values have not been normalized yet
        if contains(thisFolder,'FW','IgnoreCase',true)
            L(count,2) = 'FW';
        end
        if contains(thisFolder,'HW','IgnoreCase',true) || contains(thisFolder,'hindwing','IgnoreCase',true)
            L(count,2) = 'HW'; 
        end 
        
       if x(g) > 0
           if x(g)<1
            if y(g)<2
                L(count,2) = 'FW'; %forewing wing heart
            elseif y>7 
                L(count,2) = 'HW'; %hindwing wing heart
            end
           end
        end

        if x(g) == -0.1
            if y(g)<5
                L(count,2) = 'FWWH'; %forewing wing heart
            elseif y(g)>6 
                L(count,2) = 'HWWH'; %hindwing wing heart
            end
        end
 
        if x(g) == -2
            if y(g)>1
                L(count,2) = 'ABD'; %left and right sections of DH on ABD
            elseif y(g)<1
                L(count,2) = 'DH'; 
            end
        end       
        
        if x(g) == -0.5
            L(count,2) = 'THX'; %first thoracic segment
        elseif x(g) == -1
            L(count,2) = 'PRO'; %pronotum
        end
        
 
        count = count+1; %update L array for location
    end
    
    X(h:j,1) = x(1,:); %span values
    Y(h:j,1) = y(1,:); %chord values
    L(h:j,1) = k; %for reference to file number
    R(h:j,15) = reg; %region locale

    clear x y %otherwise code pads them onto the end
   
    h = h+m; %update h
end

%hard code in these values
L(441,2) = 'FW'; L(670,2) = 'FW'; L(671,2) = 'FW'; L(674,2) = 'FW';
L(741,2) = 'FW'; L(742,2) = 'FW';

%%
%Hard code in areas and average vein diameter (*0.5 to make into radius)
%These areas were measured in FIJI based on the high-resolution images of
%FW/HW

FW_area = 355.798*(1*10^(-6)); %FW surface area in meters squared
HW_area = 798.58*(1*10^(-6)); %HW surface area in meters squared

FW_LE_rad = 0.5*0.12*(1*10^(-3)); %FW vein radius mm to meters
FW_MEM_rad = 0.5*0.13*(1*10^(-3));
FW_TIP_rad = 0.5*0.06*(1*10^(-3));
FW_LATT_rad = 0.5*0.08*(1*10^(-3));
FW_TE_rad = 0.5*0.10*(1*10^(-3));

FW_LE_area = 63.59*(1*10^(-6));
FW_MEM_area = 9.44*(1*10^(-6));
FW_TIP_area = 47.55*(1*10^(-6));
FW_LATT_area = 149.72*(1*10^(-6));
FW_TE_area = 85.50*(1*10^(-6));

HW_LE_rad = 0.5*0.11*(1*10^(-3)); %HW vein radius mm to meters
HW_MEM_rad = 0.5*0.15*(1*10^(-3));
HW_TIP_rad = 0.5*0.04*(1*10^(-3));
HW_LATT_rad = 0.5*0.06*(1*10^(-3));
HW_TE_rad = 0.5*0.07*(1*10^(-3));

HW_LE_area = 18.75*(1*10^(-6)); %HW region surface area mm squared to meters squared
HW_MEM_area = 17.07*(1*10^(-6));
HW_TIP_area = 74.37*(1*10^(-6));
HW_LATT_area = 282.01*(1*10^(-6));
HW_TE_area = 406.48*(1*10^(-6));

wingRadii = [FW_LE_rad,FW_MEM_rad,FW_TIP_rad,FW_LATT_rad,FW_TE_rad,HW_LE_rad,HW_MEM_rad,...
    HW_TIP_rad, HW_LATT_rad, HW_TE_rad];

D_ox = 0.000018*(1*10^(-4)); % units cm2 diffusion coefficient of oxygen in water at 20.6 C (Han and Bartels 1996)

%% REORGANIZE DATA INTO MAIN PLOTTING VARS -- R, REGIONS, LOCATION, SPAN

if totalNumberOfFiles >= 1
   count = 1;
   j = 0;
   h = 1;

   for k = 1 :totalNumberOfFiles
		% Go through all those files.
		thisFolder = allFileInfo(k).folder;
		thisBaseFileName = allFileInfo(k).name;
		fullFileName = fullfile(thisFolder, thisBaseFileName);
 		fprintf('     Processing file %d of %d : "%s".\n', k, totalNumberOfFiles, fullFileName);
        
        fileParams = dir(fullfile(thisFolder, '*particle_params.csv')); %computed
        digParams = dir(fullfile(thisFolder, '*xypts.csv')); %read in digitized files
        baseParams = fullfile(thisFolder,fileParams.name); %particle analysis     
        
        if ~isempty([fileParams.name])
            C = dlmread(baseParams); %Parameter files read in 
            
            %digParams = dir(fullfile(thisFolder, '*xypts.csv')); %read in digitized files
            digParams = fullfile(thisFolder,digParams.name);

            %L(h,2) -- aka -- var is going to pass the correct region
            [partMean,numDig,partMax,partMed,numPks,cum_dist, alpha_1,re_1] = calcParticles(digParams, C(1,14), C(1,10)); %col14: physVoxSize, col10: timeInt
                
            if C(1,2) ~= 0 %if smoothed velocity data IS NOT ZERO
                r = find(C(:,2) == 0); %find indx when vel = 0 aka non-existent
                if isempty(r) 
                    vel = C(:,2); %no zeros
                    pkDist = C(:,7); %no zeros
                else
                    vel = C(1:r(1)-1,2); %index to get nonzeros
                    pkDist = C(1:r(1)-1,7); %index to get nonzeros
                end
                      
                [m,~] = size(vel); %m = number of points
                j = j+m;  %update j  
                
                R(h:j,1) = k;
                R(h:j,2) = vel; %max velocities
               % R(h:j,3) = xVel(:,1); %currently this is  mean values
                %R(h:j,4) = yVel(:,1);
                R(h:j,5) = alpha_1; %this is the first part of the Womersley equation

                R(h:j,7) = partMean; %Mean
                R(h:j,8) = partMax; %MAX
                R(h:j,9) = partMed; %Median
                R(h:j,10) = re_1; %reynolds number part 1 (everything but Length)
                R(h:j,11) = cum_dist; 
                R(h:j,12) = numDig; %number of frames digitized per particle
                R(h:j,13) = numPks; %number of pks
                R(h:j,14) =  C(1,10); %time int (sec) per slice
                %R15 = regiondata -- REDO THIS AT SOME POINT

              
                %Figure out left or right wing (done by eye previously)
                %Also differentiate between folded wings and assign L/R
                if contains(thisFolder, '2017-05-24')
                    R(h:j,6) = 0; %left wings
                    if contains(thisFolder, 'folded')
                        R(h:j,6) = 1; %right wings
                    end
                elseif contains(thisFolder, '2017-07-25')
                    R(h:j,6) = 0; %left wings  
                elseif contains(thisFolder, '2017-07-27')
                    R(h:j,6) = 0; %left wings 
                elseif contains(thisFolder, '2017-08-30')
                    if contains(thisFolder,'_right_')
                        R(h:j,6) = 1; %right wings and folded
                    else
                        R(h:j,6) = 0; %left wings
                    end
                elseif contains(thisFolder, '2017-10-26')
                    R(h:j,6) = 1; %right wings
                elseif contains(thisFolder, '2017-12-14')
                    if contains(thisFolder,'folded')
                        R(h:j,6) = 1; %R wings
                    else 
                        R(h:j,6) = 0; %L wings
                    end
                elseif contains(thisFolder, '2017-12-15')
                    if contains(thisFolder,'folded')
                        R(h:j,6) = 0; %left wings
                    else
                        R(h:j,6) = 1; %right wings
                    end
                elseif contains(thisFolder, '2017-12-18')
                    if contains(thisFolder,'_right_')
                        R(h:j,6) = 1; %right wings
                    elseif contains(thisFolder,'_left_')
                        R(h:j,6) = 0; %left wings
                    end
                end
                h = h+m; %update h
            end           
        end        
		[~, baseNameNoExt, ~] = fileparts(thisBaseFileName);
		clear index
        count = count + 1;
    end
    
k = L(:,2);
j = contains(k,'ABD','IgnoreCase',true); R(j,6) = -2.5;   
j = contains(k,'PRO','IgnoreCase',true); R(j,6) = -1;
j = contains(k,'THX','IgnoreCase',true); R(j,6) = -0.5;  
j = contains(k,'DH','IgnoreCase',true);  R(j,6) = -2;
j = contains(k,'FWWH','IgnoreCase',true); R(j,6) = -0.10;
j = contains(k,'HWWH','IgnoreCase',true); R(j,6) = -0.15;
else
	fprintf('     Folder %s has no files in it.\n', thisFolder);
end

fprintf('\nDone looking in all %d folders!\nFound %d files in the %d folders.\n', numberOfFolders, totalNumberOfFiles, numberOfFolders);
%% COLOR for all metrics - vel, frequency, path sinu, net flow
color_path = '/Users/marysalcedo/Documents/MATLAB/Wings_Matlab_Genscode/';
addpath('/Users/marysalcedo/Documents/MATLAB/Wings_Matlab_Genscode/cbrewer');
%80 was chosen based on vector size of 788 particles
%By setting the last number, you create a 3 column color vector with that many rows
colorF = cbrewer('seq', 'BuPu',80); %colormap length of A
colorD = cbrewer('seq', 'BuGn',80); %colormap length of A
colorM = cbrewer('seq', 'PuRd',80); %colormap length of A
colorB = cbrewer('div', 'RdYlBu',80); %colormap length of A (Velocity)
colorFreq = cbrewer('seq', 'OrRd',80); %colormap for freq

%Create a more targeted wing color spectrum. Smaller range will show
%differences between individual particles better than the size 80
colorBsmall = cbrewer('div', 'RdYlBu',59); %colormap length of A (Velocity)
flipcolorB = flip(colorBsmall); %on 18Dec2019, flipping colors so it runs cold to hot

colorDsmall = cbrewer('seq', 'BuGn',59); %colormap length of path sinu
flipcolorD = colorDsmall; %on 18Dec2019, NO FLIP, but leaving name flipping colors so it runs cold to hot

colorMsmall = cbrewer('seq', 'PuRd',59); %colormap length of net flow
flipcolorM = colorMsmall; %on 18Dec2019, NO FLIP but leaving name flipping colors so it runs cold to hot

colorF2 = cbrewer('seq', 'BuPu',10); %colormap length of A
colorD2 = cbrewer('seq', 'BuGn',10); %colormap length of A
colorB2 = cbrewer('div', 'RdYlBu',10); %colormap length of A
colorB2 = flip(colorB2);
colorM2 = cbrewer('seq', 'PuRd',10); %colormap length of A
colorFreq2 = cbrewer('seq', 'OrRd',10); %colormap for freq

%Made on coolors.co/palettes/trending
colorKeyLimeBlue = [217,237,146;...
             181, 228, 140;...
             153, 217, 140;...
             118, 200, 147;...
             82, 182, 154;...
             52, 160, 164;...
             22, 138, 173;...
             26, 117, 159;...
             30, 96, 145;...
             24, 78, 119;]/256;
         
         colorKeyLimeBlue8 = [217,237,146;...
             153, 217, 140;...
             118, 200, 147;...
             82, 182, 154;...
             52, 160, 164;...
             22, 138, 173;...
             26, 117, 159;...
             30, 96, 145]/256;
         
          colorKeyLimeBlue6 = [217,237,146;...
             153, 217, 140;...
             118, 200, 147;...
             82, 182, 154;...
             22, 138, 173;...
             26, 117, 159]/256;
         




%% Separate FW, HW, heart, thorax and abd
speedVar = 7; %7 = mean, 8 = max, 9 = min -- change this to get at diff velocities
distVar = 11;

%Volumetric flow rate -- MEANDERING
wingME = R(:,11)./(R(:,12).*R(:,14)); %volume rate = cumulative distance/(frames*timeInt)

%Entire wing
wingMEVAR = NaN(length(R),1); %wing meandering
%wingDist = NaN(length(R),4); %within main block of code
wingPulse = NaN(length(R),7); % for entire FW or HW -- within main block of code
wingWo = NaN(length(R),1); %why 7??
wingRE =  NaN(length(R),1); %wing reynolds number (will reference column R10)
wingCumDist = NaN(length(R),1); 
wingP = NaN(length(R),1); 
wingMed = NaN(length(R),1); 
wingPe = NaN(length(R),1); %peclet number

%By regions
regions = NaN(length(R),28); %predefine
regionsMedian = NaN(length(R),28); %predefine for median velocities R col9
wingReg = NaN(length(R),10); %REG SPEED (within main block of code)
wingRegMax = NaN(length(R),10); 
wingRegFreq = NaN(length(R),10); %REG FREQ
wingRegMe = NaN(length(R),10); %REG MEANDERING
wingRegDIST = NaN(length(R),10);
wingRegWo = NaN(length(R),10); %predefine the Womersley number (part 2)
wingRegRE =  NaN(length(R),10); %wing reynolds number by region (will reference column R10)
wingRegMed =  NaN(length(R),10); 
wingRegPe  = NaN(length(R),10); %peclet number


wingMax = NaN(length(R),1); %Sort max speeds for just fw/hw

%wingDist2 = NaN(length(R),1); %just use to create sorted list


%Assign variables
%Matrix R, column 6: 1s and 0s indicate FW or HW and regions in the body have negative values (Rcol6 defined
%above). If data is not grabbed, it is replaced by an NaN.
%O and 1 indicates Left and Right
for i = 1:length(R)
    if R(i,6) == 1 %Grab all FW and HW data
        wingP(i,1) = R(i,13)/(R(i,12)*R(i,14)); %freq = num peaks/(time)
        wingMax(i,1) = R(i,8); %max speeds for FW/HW
        %wingDist2(i,1) = R(i,10); % start to end distance
        wingCumDist(i,1) = R(i,11); %cumulative dist
        wingMEVAR(i,1) = wingME(i,1); %wing meandering
        wingWo(i,1) = R(i,5); %col5 contains Womersley part 1
        wingRE(i,1) = R(i,10); %col10 contains Reynolds part 1: rho*V/mu
        wingMed(i,1) = R(i,9); %col9 contains median velocities
    end
    %AND grab values at 0
    if R(i,6) == 0 %Grab all FW and HW data
        wingP(i,1) = R(i,13)/(R(i,12)*R(i,14));
        wingMax(i,1) = R(i,8); %max speeds for FW/HW
        %wingDist2(i,1) = R(i,10); % start to end distance
        wingCumDist(i,1) = R(i,11); %cumulative dist
        wingMEVAR(i,1) = wingME(i,1); %wing meandering
        wingWo(i,1) = R(i,5); %col5 contains Womersley part 1
        wingRE(i,1) = R(i,10); %col10 contains Reynolds part 1: rho*V/mu
        wingMed(i,1) = R(i,9); %col9 contains median velocities

    end
end

%Cumulative/total distance
% OLD: distTRAV = wingCumDist./wingDist2; %PARAMETER FOR wing-only DIST
distTRAV = wingCumDist; %new total distance traveled
distTRAV(~isfinite(distTRAV)) = 0; %replace inf with 0's
sortTRAV = sort(distTRAV); 
colorTRAVH = NaN(length(R),3); colorTRAVF = NaN(length(R),3);


%Peak Speeds
sortAllPeakSpeeds = sort(R(:,8)); %sort all peak speeds
sortWingMax = sort(wingMax); %629 -THIS IS THE COLOR COMPARISON FOR FW/HW

%Volume Rate - Cumulate distance over time
sortWingME = sort(wingMEVAR); %this is for color comp
wingMEFW = NaN(length(R),1);           wingMEHW = NaN(length(R),1);
colorMEFW = NaN(length(R),3);          colorMEHW = NaN(length(R),3);

colorspeedF = NaN(length(R),3);        colorspeedH = NaN(length(R),3);
colorspeedABD = NaN(length(R),3);      colorspeedDH = NaN(length(R),3);
colorspeedPRO = NaN(length(R),3);      colorspeedTHX = NaN(length(R),3);
colorspeedFWWH = NaN(length(R),3);     colorspeedHWWH = NaN(length(R),3);
colorspeedmeanF = NaN(length(R),3);    colorspeedmeanH = NaN(length(R),3);

%Particle Oscillation  
sortP = sort(wingP); %WING FREQ
colorFreqPFW =  NaN(length(wingP),3); colorFreqPHW =  NaN(length(wingP),3);

%Womersley
%don't worry about color yet

%For the whole length of the data
count = 0;
for k = 1:length(X)
    if R(k,12)>25 %cut out low frame number digitizations
        
        %This is for color sepctrum
        val = R(k,speedVar); %speed var for all other regions
        [~,a] = find(val <= sortAllPeakSpeeds(:,1)); %find index of speedvalue
        if ~isnan(a)
            indx = round(a(1)/10); %since high #s, divide and round up
            if indx == 0 
                indx = a(1);
            end
        end
        
        %FW/HW color index
        val2 = R(k,speedVar);
        a = find(val2 == sortWingMax(:,1));
        if ~isnan(a)
            count = count+1;
            indx2 = round(a(1)/10); %since high #s, divide and round up
            if indx2 == 0 
                indx2 = a(1);
            end
        end
        val3 = distTRAV(k,1); % color index CUMULATIVE DIST
        a = find(val3 <= sortTRAV(:,1));
        if ~isnan(a)
            indx3 = round(a(1)/10); %since high #s, divide and round up
            if indx3 == 0 
                indx3 = a(1);
            end
        end
        %colorindex
        val4 = wingMEVAR(k,1); %MEANDERING ONLY WINGS- NET VOLUME
        a = find(val4 <= sortWingME(:,1));
        if ~isnan(a)
            indx4 = round(a(1)/10); %since high #s, divide and round up
            if indx4 == 0 
                indx4 = a(1);
            end
        end
        
        %color index
        val5 = wingP(k,1); %wing freq -- this is all for the color palette
        if val5 == 0 
            indx5 = 1;
        else
            a = find(val5 <= sortP(:,1));
            if ~isnan(a)
                indx5 = round(a(1)/10); %since high #s, divide and round up
                if indx5 == 0 
                    indx5 = a(1);
                end
            end
        end
        

        %separate the regions
        %L is a vector that contains BIG regions: FW, FWWH, HW, HWWH, ABD,
        %THX, PRO, and DH. This is more reliable than probing column 15 of
        %R which also has regions of 0 or NaN in it.
        if contains(L(k,2),'FWWH','IgnoreCase',true) %if FW, grab data
            regions(k,1) = R(k,speedVar); %take data point from R
            regionsMedian(k,1) = R(k,9); %col9 is median velocities
            regions(k,2) = X(k); %take span loc from above
            regions(k,3) = Y(k); %take chord loc from above
          
            wingPulse(k,3) =  (R(k,13))/(R(k,12)*R(k,14)); %pks*numFrames/(time_int*numFrames)

            %for each value, grab approp. color index from set of 80 colors
            colorspeedFWWH(k,:) = colorB(indx,:); %whole row
        elseif contains(L(k,2),'HWWH','IgnoreCase',true)
            regions(k,4) = R(k,speedVar); %take data point from R
            regionsMedian(k,4) = R(k,9); %col9 is median velocities
            regions(k,5) = X(k); %take span loc from R
            regions(k,6) = Y(k); %take chord loc from R
            wingPulse(k,4) =  (R(k,13))/(R(k,12)*R(k,14)); %pks*numFrames/(time_int*numFrames)

            colorspeedHWWH(k,:) = colorB(indx,:);
        elseif contains(L(k,2),'ABD','IgnoreCase',true) 
            regions(k,7) = R(k,speedVar); %take data point from R
            regionsMedian(k,7) = R(k,9); %col9 is median velocities
            regions(k,8) = X(k); %take span loc from R
            regions(k,9) = Y(k); %take chord loc from R
            wingPulse(k,7) =  (R(k,13))/(R(k,12)*R(k,14)); %pks*numFrames/(time_int*numFrames)
            colorspeedABD(k,:) = colorB(indx,:);
        elseif contains(L(k,2),'DH','IgnoreCase',true)
            regions(k,10) = R(k,speedVar); %take data point from R
            regionsMedian(k,10) = R(k,9); %col9 is median velocities
            regions(k,11) = X(k); %take span loc from R
            regions(k,12) = Y(k); %take chord loc from R
            
            wingPulse(k,5) =  (R(k,13))/(R(k,12)*R(k,14)); %pks*numFrames/(time_int*numFrames)

            colorspeedDH(k,:) = colorB(indx,:);
        elseif contains(L(k,2),'THX','IgnoreCase',true) 
            regions(k,13) = R(k,speedVar); %take data point from R
            regionsMedian(k,13) = R(k,9); %col9 is median velocities
            regions(k,14) = X(k); %take span loc from R
            regions(k,15) = Y(k); %take chord loc from R
            
            wingPulse(k,6) =  (R(k,13))/(R(k,12)*R(k,14)); %pks*numFrames/(time_int*numFrames)
            colorspeedTHX(k,:) = colorB(indx,:);
        elseif contains(L(k,2),'PRO','IgnoreCase',true) 
            regions(k,16) = R(k,speedVar); %take data point from R
            regionsMedian(k,16) = R(k,9); %col9 is median velocities
            regions(k,17) = X(k); %take span loc from R
            regions(k,18) = Y(k); %take chord loc from R
            colorspeedPRO(k,:) = colorB(indx,:);
        elseif contains(L(k,2),'FW','IgnoreCase',true) 
            regions(k,19) = R(k,speedVar); %take data point from R
            regionsMedian(k,19) = R(k,9); %col9 is median velocities
            regions(k,20) = R(k,3);
            regions(k,21) = R(k,4);
            regions(k,22) = X(k); %take span loc 
            regions(k,23) = Y(k); %take chord loc 
            
            wingDist(k,1) = distTRAV(k,1); %path sinu
            wingMEFW(k,1) = wingMEVAR(k,1); %net flow
            wingPulse(k,1) = wingP(k,1); %pks*numFrames/(time_int*numFrames)
            %wingWo(k,1) = ??; Setting up Womersley for whole wing

            
            if R(k,15)== 1
                wingReg(k,1) = R(k,speedVar);
                wingRegMax(k,1) = wingMax(k,1);
                wingRegFreq(k,1) = wingP(k,1);
                wingRegMe(k,1) = wingMEVAR(k,1);
                wingRegDIST(k,1) = distTRAV(k,1);
                wingRegWo(k,1) = wingWo(k,1)*FW_LE_rad;
                wingRegRE(k,1) = wingRE(k,1)*FW_LE_rad;
                wingRegMed(k,1) = wingMed(k,1);
                wingRegPe(k,1) = (wingMed(k,1)*(1*10^(-6))*FW_LE_rad)/D_ox;
                
            elseif R(k,15)== 2
                wingReg(k,2) = R(k,speedVar);
                wingRegMax(k,2) = wingMax(k,1);
                wingRegFreq(k,2) = wingP(k,1);
                wingRegMe(k,2) = wingMEVAR(k,1);
                wingRegDIST(k,2) = distTRAV(k,1);
                wingRegWo(k,2) =  wingWo(k,1)*FW_MEM_rad;
                wingRegRE(k,2) =  wingRE(k,1)*FW_MEM_rad;
                wingRegMed(k,2) = wingMed(k,1);
                wingRegPe(k,2) = (wingMed(k,1)*(1*10^(-6))*FW_MEM_rad)/D_ox;

            elseif R(k,15)== 3
                wingReg(k,3) = R(k,speedVar);
                wingRegMax(k,3) = wingMax(k,1);
                wingRegFreq(k,3) = wingP(k,1);
                wingRegMe(k,3) = wingMEVAR(k,1);
                wingRegDIST(k,3) = distTRAV(k,1);
                wingRegWo(k,3) =  wingWo(k,1)*FW_TIP_rad;
                wingRegRE(k,3) =  wingRE(k,1)*FW_TIP_rad;
                wingRegMed(k,3) = wingMed(k,1);
                wingRegPe(k,3) = (wingMed(k,1)*(1*10^(-6))*FW_TIP_rad)/D_ox;

            elseif R(k,15)== 4
                wingReg(k,4) = R(k,speedVar);
                wingRegMax(k,4) = wingMax(k,1);
                wingRegFreq(k,4) = wingP(k,1);
                wingRegMe(k,4) = wingMEVAR(k,1);
                wingRegDIST(k,4) = distTRAV(k,1);
                wingRegWo(k,4) =  wingWo(k,1)*FW_LATT_rad;
                wingRegRE(k,4) =  wingRE(k,1)*FW_LATT_rad;
                wingRegMed(k,4) = wingMed(k,1);
                wingRegPe(k,4) = (wingMed(k,1)*(1*10^(-6))*FW_LATT_rad)/D_ox;
                
            elseif R(k,15)== 5
                wingReg(k,5) = R(k,speedVar);
                wingRegMax(k,5) = wingMax(k,1);
                wingRegFreq(k,5) = wingP(k,1);
                wingRegMe(k,5) = wingMEVAR(k,1);
                wingRegDIST(k,5) = distTRAV(k,1);
                wingRegWo(k,5) =  wingWo(k,1)*FW_TE_rad;
                wingRegRE(k,5) =  wingWo(k,1)*FW_TE_rad;
                wingRegMed(k,5) = wingMed(k,1);
                wingRegPe(k,5) = (wingMed(k,1)*(1*10^(-6))*FW_TE_rad)/D_ox;
            end
            
            colorFreqPFW(k,:) = colorF(indx5,:); %All I care about is color indx
            colorspeedF(k,:) = colorB(indx2,:); %based on all peaks speeds
            colorTRAVF(k,:) = colorD(indx3,:); %based on all peaks speeds
            colorMEFW(k,:) = colorM(indx4,:);

            %colorspeedmeanF(k,:) = colorB(indx2(1),:); %based on ONLY FW/HW peak speeds  
            
        elseif contains(L(k,2),'HW','IgnoreCase',true) 
            regions(k,24) = R(k,speedVar); %take data point from R
            regionsMedian(k,24) = R(k,9); %col9 is median velocities
            regions(k,25) = R(k,3); %x vels
            regions(k,26) = R(k,4); %y vels
            regions(k,27) = X(k); %take span loc from R
            regions(k,28) = Y(k); %take chord loc from R

            wingDist(k,2) = distTRAV(k,1); %path sinu
            wingMEHW(k,1) = wingMEVAR(k,1);
            wingPulse(k,2) = wingP(k,1); %pks*numFrames/(time_int*numFrames)
            
            if R(k,15)== 1
                wingReg(k,6) = R(k,speedVar);
                wingRegMax(k,6) = wingMax(k,1);
                wingRegFreq(k,6) = wingP(k,1);
                wingRegMe(k,6) = wingMEVAR(k,1);
                wingRegDIST(k,6) = distTRAV(k,1);
                wingRegWo(k,6) =  wingWo(k,1)*HW_LE_rad;
                wingRegRE(k,6) =  wingRE(k,1)*HW_LE_rad;
                wingRegMed(k,6) = wingMed(k,1);
                wingRegPe(k,6) = (wingMed(k,1)*(1*10^(-6))*HW_LE_rad)/D_ox;
                
            elseif R(k,15)== 2
                wingReg(k,7) = R(k,speedVar); 
                wingRegMax(k,7) = wingMax(k,1);
                wingRegFreq(k,7) = wingP(k,1);
                wingRegMe(k,7) = wingMEVAR(k,1);
                wingRegDIST(k,7) = distTRAV(k,1);
                wingRegWo(k,7) =  wingWo(k,1)*HW_MEM_rad;
                wingRegRE(k,7) =  wingRE(k,1)*HW_MEM_rad;
                wingRegMed(k,7) = wingMed(k,1);
                wingRegPe(k,7) = (wingMed(k,1)*(1*10^(-6))*HW_MEM_rad)/D_ox;

                
            elseif R(k,15)== 3
                wingReg(k,8) = R(k,speedVar);
                wingRegMax(k,8) = wingMax(k,1);
                wingRegFreq(k,8) = wingP(k,1);
                wingRegMe(k,8) = wingME(k,1);
                wingRegDIST(k,8) = distTRAV(k,1);
                wingRegWo(k,8) =  wingWo(k,1)*HW_TIP_rad;
                wingRegRE(k,8) =  wingRE(k,1)*HW_TIP_rad;
                wingRegMed(k,8) = wingMed(k,1);
                wingRegPe(k,8) = (wingMed(k,1)*(1*10^(-6))*HW_TIP_rad)/D_ox;
                
            elseif R(k,15)== 4
                wingReg(k,9) = R(k,speedVar);
                wingRegMax(k,9) = wingMax(k,1);
                wingRegFreq(k,9) = wingP(k,1);
                wingRegMe(k,9) = wingMEVAR(k,1);
                wingRegDIST(k,9) = distTRAV(k,1);
                wingRegWo(k,9) =  wingWo(k,1)*HW_LATT_rad;
                wingRegRE(k,9) =  wingRE(k,1)*HW_LATT_rad;
                wingRegMed(k,9) = wingMed(k,1);
                wingRegPe(k,9) = (wingMed(k,1)*(1*10^(-6))*HW_LATT_rad)/D_ox;
                
            elseif R(k,15)== 5
                wingReg(k,10) = R(k,speedVar);
                wingRegMax(k,10) = wingMax(k,1);
                wingRegFreq(k,10) = wingP(k,1);
                wingRegMe(k,10) = wingMEVAR(k,1);
                wingRegDIST(k,10) = distTRAV(k,1);
                wingRegWo(k,10) =  wingWo(k,1)*HW_TE_rad;
                wingRegRE(k,10) =  wingRE(k,1)*HW_TE_rad;
                wingRegMed(k,10) = wingMed(k,1);
                wingRegPe(k,10) = (wingMed(k,1)*(1*10^(-6))*HW_TE_rad)/D_ox;
                
            end
            
            colorspeedH(k,:) = colorB(indx2,:);%based on all wing peaks speeds
            colorTRAVH(k,:) = colorD(indx3,:);%path sinu
            colorMEHW(k,:) = colorM(indx4,:);%net flow
            colorFreqPHW(k,:) = colorF(indx5,:);%freq

            %colorspeedmeanH(k,:) = colorB(indx2(1),:);%based on ONLY FW/HW peak speeds   
        end  
    end
end
%% 18DEC2019 - new color palette based on 333 particles for FW, 252 for HW

% flipcolorB is set to have 34 values

fwEDIT = regions(:,19); 
%fwEDIT = wingMEFW(:,1); % for FW net flow
%fwEDIT = wingDist(:,1); %for path sinu
fwEDIT_out = fwEDIT(~isnan(fwEDIT)); %remove NaNs from forewing vel vector


%HW has 252 particles 
hwEDIT = regions(:,24); %vel
%hwEDIT = wingMEHW(:,1); %for net flow
%hwEDIT = wingDist(:,2); %for path sinu
hwEDIT_out = hwEDIT(~isnan(hwEDIT)); %remove NaNs from hindwing vel vector


combo = [fwEDIT_out;hwEDIT_out]; %combine FW/HW vectors into one (end-on-end)
sortcombo = sort(combo); %sorted list of all velocity values in FW/HW

fwEDIT_X = regions(:,22); 
fwEDIT_Y = regions(:,23); 
fwEDIT_outX = fwEDIT_X(~isnan(fwEDIT_X)); %remove NaNs from forewing Xs
fwEDIT_outY = fwEDIT_Y(~isnan(fwEDIT_Y)); %remove NaNs from forewing Ys

fwEDIT_outSORT = sort(fwEDIT_out); %sort FW velocity values
colorspeedF_EDIT = NaN(length(fwEDIT_out),3); 

for i = 1:length(fwEDIT_out)
    indx(i) = find(fwEDIT_out(i) == sortcombo); %vector of indices -- compare to ALL vels
    
    indx2(i) = round(indx(i)/10); %since high #s, divide and round up
    if indx2(i) == 0 
        indx2(i) = indx(1,i); %means indx(i) is under 10 
    end
    %change to D for paht sinu
    %change to M for net flow
    %change to B for vel
    colorspeedF_EDIT(i,:) = flipcolorB(indx2(i),:); 
end

%THIS ISN'T PRETTY -- BUT I JUST WANT TO PLOT, SO COPY/PASTE
hwEDIT_X = regions(:,27); 
hwEDIT_Y = regions(:,28); 
hwEDIT_outX = hwEDIT_X(~isnan(hwEDIT_X)); %remove NaNs from forewing Xs
hwEDIT_outY = hwEDIT_Y(~isnan(hwEDIT_Y)); %remove NaNs from forewing Ys

hwEDIT_outSORT = sort(hwEDIT_out); %sort velocity values
colorspeedH_EDIT = NaN(length(hwEDIT_out),3); 

%indices had to either be cleared or called something diff
for i = 1:length(hwEDIT_out)
    indx3(i) = find(hwEDIT_out(i) == sortcombo); %vector of indices  -- compare to ALL vels
    
    indx4(i) = round(indx3(1,i)/10); %since high #s, divide and round up
    if indx4(i) == 0 
        indx4(i) = indx3(1,i); %means indx(i) is under 10 
    end
    colorspeedH_EDIT(i,:) = flipcolorB(indx4(i),:);
end



%%  REGION ANALYSIS - PULSE FREQUENCY
[~,n] = size(wingReg);
figure
for i = 1:n %go through fw and hw
    freqReg(i,1) = median(wingRegFreq(:,i),'omitnan')
end
sortFreqReg = sort(freqReg(:,1)); %sort wing speeds
x = 1:10; %each of the wing regions
r = sortFreqReg; 

for i = 1:n
    indx = find(freqReg(i,1) == r(:,1)); %find when index of the speed matches the sorted speed
    scatter(x(i),freqReg(i,1),900, colorKeyLimeBlue(indx(1),:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
    hold on
end
colormap(colorKeyLimeBlue)
c = colorbar('Ticks', [0,1], 'TickLabels', {min(r(:,1)), max(r(:,1))}, 'FontSize', 14);
c.Label.String = 'Mean frequency';
c.Label.FontSize = 14;
title('Color palette of frequency')
hold off
%%  REGION ANALYSIS - AVERAGE VEIN DIAMETERS/RADII
[~,n] = size(wingReg);
figure

sortRadReg = sort(wingRadii(1,:)); %sort wing lengths
x = 1:10; %each of the wing regions
r = sortRadReg; 

for i = 1:n
    indx = find(wingRadii(1,i) == r(1,:)); %find when index fof the length matches the sorted length
    scatter(x(i),wingRadii(1,i),900, colorKeyLimeBlue(indx(1),:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
    hold on
end
colormap(colorKeyLimeBlue)
c = colorbar('Ticks', [0,1], 'TickLabels', {min(r(1,:)), max(r(1,:))}, 'FontSize', 14);
c.Label.String = 'Average Vein diameter (m)';
c.Label.FontSize = 14;
title('Color palette of avg. vein diameter')
hold off
%%  REGION ANALYSIS - REYNOLDS
[~,n] = size(wingReg); %Go through 10 columns, all the same

figure
for i = 1:n %go through fw and hw
    REReg(i,1) = median(wingRegRE(:,i),'omitnan')
end
sortREReg = sort(REReg(:,1)); %sort wing speeds
x = 1:10; %each of the wing regions
r = sortREReg 

for i = 1:n
    indx = find(REReg(i,1) == r(:,1)); %find when index of the speed matches the sorted speed
    scatter(x(i),REReg(i,1),900, colorKeyLimeBlue(indx(1),:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
    hold on
end
colormap(colorKeyLimeBlue)
c = colorbar('Ticks', [0,1], 'TickLabels', {min(r(:,1)), max(r(:,1))}, 'FontSize', 14);
c.Label.String = 'REYNOLDS NUMBER';
c.Label.FontSize = 14;
title('Color palette of MEAN RE')
hold off
%% REGION ANALYSIS - MEDIAN VELOCITIES
[~,n] = size(wingReg); %Go through 10 columns, all the same
figure
for i = 1:n %go through fw and hw
    medReg(i,1) = mean(wingRegMed(:,i),'omitnan')
end
sortMedReg = sort(medReg(:,1)); %sort wing speeds
x = 1:10; %each of the wing regions
r = sortMedReg; 

for i = 1:n
    indx = find(medReg(i,1) == r(:,1)); %find when index of the speed matches the sorted speed
    scatter(x(i),medReg(i,1),900, colorKeyLimeBlue(indx(1),:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
    hold on
end
colormap(colorKeyLimeBlue)
c = colorbar('Ticks', [0,1], 'TickLabels', {min(r(:,1)), max(r(:,1))}, 'FontSize', 14);
c.Label.String = 'Median Velocities';
c.Label.FontSize = 14;
title('Median Velocities')
hold off
%% REGION ANALYSIS - PECLET NUMBER
[~,n] = size(wingReg); %Go through 10 columns, all the same
figure
for i = 1:n %go through fw and hw
    medPe(i,1) = mean(wingRegPe(:,i),'omitnan')
end
sortMedPe = sort(medPe(:,1)); %sort wing speeds
x = 1:10; %each of the wing regions
r = sortMedPe; 

for i = 1:n
    indx = find(medPe(i,1) == r(:,1)); %find when index of the speed matches the sorted speed
    scatter(x(i),medPe(i,1),900, colorKeyLimeBlue(indx(1),:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
    hold on
end
colormap(colorKeyLimeBlue)
c = colorbar('Ticks', [0,1], 'TickLabels', {min(r(:,1)), max(r(:,1))}, 'FontSize', 14);
c.Label.String = 'Peclet';
c.Label.FontSize = 14;
title('Peclet')
hold off
%% REGION ANALYSIS - MAX VELOCITIES
[~,n] = size(wingReg);
figure
for i = 1:n %go through fw and hw
    speedReg(i,1) = mean(wingReg(:,i),'omitnan')
end
sortspeedReg = sort(speedReg(:,1)); %sort wing speeds
x = 1:10; %each of the wing regions
r = sortspeedReg; 

for i = 1:n
    indx = find(speedReg(i,1) == r(:,1)); %find when index of the speed matches the sorted speed
    scatter(x(i),speedReg(i,1),900, colorKeyLimeBlue(indx(1),:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
    hold on
end
colormap(colorKeyLimeBlue)
c = colorbar('Ticks', [0,1], 'TickLabels', {min(r(:,1)), max(r(:,1))}, 'FontSize', 14);
c.Label.String = 'Mean Speed (micon/sec)';
c.Label.FontSize = 14;
title('Color palette of mean, regional speed')
hold off
%% REGION ANALYSIS - MAX VELOCITIES TEST CASE VERSION 2
% It performs same as wingReg. I was just checking
[~,n] = size(wingReg);
figure
for i = 1:n %go through fw and hw
    speedRegMax(i,1) = mean(wingRegMax(:,i),'omitnan')
end
sortspeedRegMax = sort(speedRegMax(:,1)); %sort wing speeds
x = 1:10; %each of the wing regions
r = sortspeedRegMax; 

for i = 1:n
    indx = find(speedRegMax(i,1) == r(:,1)); %find when index of the speed matches the sorted speed
    scatter(x(i),speedRegMax(i,1),900, colorKeyLimeBlue(indx,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
    hold on
end
colormap(colorKeyLimeBlue)
c = colorbar('Ticks', [0,1], 'TickLabels', {min(r(:,1)), max(r(:,1))}, 'FontSize', 14);
c.Label.String = 'Max Velocities...test case';
c.Label.FontSize = 14;
title('Color palette of max velocities')
hold off
%% REGION ANALYSIS - WOMERSLEY
figure
[~,n] = size(wingRegWo);
for i = 1:n %go through fw and hw
    WoReg(i,1) = median(wingRegWo(:,i),'omitnan')
end
sortWoReg = sort(WoReg(:,1)); %sort wing Wo
x = 1:10;
r = sortWoReg;

for i = 1:n
    indx = find(WoReg(i,1) == sortWoReg(:,1));
    scatter(x(i),WoReg(i,1),900, colorKeyLimeBlue(indx(1),:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
    hold on
end
colormap(colorKeyLimeBlue)
c = colorbar('Ticks', [0,1], 'TickLabels', {min(r(:,1)), max(r(:,1))}, 'FontSize', 14);
c.Label.String = 'Womersley';
c.Label.FontSize = 14;
title('Color palette of regional Womersley')
hold off
%% REGION ANALYSIS - INSTANTANEOUS SPEED
figure
for i = 1:n %go through fw and hw
    MEReg(i,1) = median(wingRegMe(:,i),'omitnan')
end
sortMEReg = sort(MEReg(:,1)); %
x = 1:10;
r = sortMEReg;

for i = 1:n
    indx = find(MEReg(i,1) == sortMEReg(:,1));
    scatter(x(i),MEReg(i,1),900, colorKeyLimeBlue(indx,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
    hold on
end
colormap(colorKeyLimeBlue)
c = colorbar('Ticks', [0,1], 'TickLabels', {min(r(:,1)), max(r(:,1))}, 'FontSize', 14);
c.Label.String = 'Velocity over summed distance (m/s)';
c.Label.FontSize = 14;
title('Color palette of median, instantaneous speed')
hold off

%% REGION ANALYSIS - CUMULATIVE DISTANCE
figure
n = 10;
for i = 1:n %go through fw and hw
    distReg(i,1) = median(wingRegDIST(:,i),'omitnan') %pick out regional means, maxes
end
sortdistReg = sort(distReg(:,1)); %sort CUMULATIVE DISTANCE
x = 1:10;
r = sortdistReg;

for i = 1:n
    indx = find(distReg(i,1) == r(:,1));
    scatter(x(i),distReg(i,1),900, colorD2(indx,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
    hold on
end
colormap(colorD2)
c = colorbar('Ticks', [0,1], 'TickLabels', {min(r(:,1)), max(r(:,1))}, 'FontSize', 14);
c.Label.String = 'Particle path sinuosity';
c.Label.FontSize = 14;
title('Color palette of mean - cumulative distance')
hold off

%% BODY METRICS

%Category order:'ABD','PRO','DH','SCUT','FW Hinge','FW','HW Hinge','HW'

%Take mean of peak speeds
bodyMaxMean = mean([regions(:,7),regions(:,16),regions(:,10),...
    regions(:,13),regions(:,1),regions(:,19),regions(:,4),regions(:,24)], 'omitnan');

%Take median of velocity max
bodyMaxMed = median([regions(:,7),regions(:,16),regions(:,10),...
    regions(:,13),regions(:,1),regions(:,19),regions(:,4),regions(:,24)], 'omitnan');


%Take mean of median velocities
bodyAvgMedian = mean([regionsMedian(:,7),regionsMedian(:,16),regionsMedian(:,10),...
    regionsMedian(:,13),regionsMedian(:,1),regionsMedian(:,19),regionsMedian(:,4),regionsMedian(:,24)], 'omitnan');


%Take median of velocity medians
bodyMedOfMed = median([regionsMedian(:,7),regionsMedian(:,16),regionsMedian(:,10),...
    regionsMedian(:,13),regionsMedian(:,1),regionsMedian(:,19),regionsMedian(:,4),regionsMedian(:,24)], 'omitnan');

%% BODY PULSE FREQUENCY METRICS

%Category order:'DH','SCUT','FW Hinge','FW','HW Hinge','HW'
bodyFreqMean = mean([wingPulse(:,5),wingPulse(:,6),wingPulse(:,3),wingPulse(:,1),...
    wingPulse(:,4),wingPulse(:,2)],'omitnan');


%%  BODY REGION ANALYSIS - BODY MAX VELOCITIES

close all;
[~,n] = size(bodyFreqMean); % There are 8 body regions
figure

n = 6; %if doing body frequency use this
var = bodyFreqMean;

sortVar = sort(var(1,:)); %sort body velocities (max)
x = 1:n; %each of the body regions
r = sortVar; 

for i = 1:n
    indx = find(var(1,i) == r(1,:)); %find when index fof the length matches the sorted length
    scatter(x(i),var(1,i),900, colorKeyLimeBlue8(indx(1),:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
    hold on
end
colormap(colorKeyLimeBlue6)
c = colorbar('Ticks', [0,1], 'TickLabels', {min(r(1,:)), max(r(1,:))}, 'FontSize', 14);
c.Label.String = 'Body Frequencies - Mean';
c.Label.FontSize = 14;
title('Color palette of body freq mean')
hold off

%% BOXPLOTS OF MEDIAN VEL AND STATS
%close all
figure
set(gcf,'color','w') %change whole figure background to white

labels = {'ABD','PRO','DH','SCUT','FW Hinge','HW Hinge','FW','HW'}; %labels for box plots

%Use of CategoricalScatter to plot the columns of my data
CategoricalScatterplot([regionsMedian(:,7),regionsMedian(:,16),regionsMedian(:,10),regionsMedian(:,13),...
    regionsMedian(:,1),regionsMedian(:,4),regionsMedian(:,19),regionsMedian(:,24)],...
    'MedianLineWidth', 6, 'MarkerSize',50, 'WhiskerLineWidth', 4,...
    'WhiskerColor', [0.3438 0.3477 0.3477],'MedianColor',[0.4492 0.0469 0.0078],...
    'BoxColor',[0.6523 0.6602 0.6719],'BoxAlpha',1,'Labels',labels)

%Appropriate labels
ylabel('Median velocity','FontSize',20)

xlabel('Wing and Body Regions','FontSize',20)
%axis([0 8 -10 80])
ax = gca; ax.FontSize = 14; %change axes numbers to larger font

%This is if I want stats...
dat = log([regionsMedian(:,7),regionsMedian(:,16),regionsMedian(:,10),...
    regionsMedian(:,13),regionsMedian(:,1),...
    regionsMedian(:,4),regionsMedian(:,19),regionsMedian(:,24)]);

[~,~,stats] = anova1(dat);
c = multcompare(stats)
dlmwrite('stats_medianV_Pvals.csv',c)
%%



statsRegions = zeros(6,22); %create emtpy array - 22 = time chunks

reg = log([regionsMedian(:,7),regionsMedian(:,16),regionsMedian(:,10),regionsMedian(:,13),...
    regionsMedian(:,1),regionsMedian(:,19),regionsMedian(:,4),regionsMedian(:,24)]);

for i = 1:2%length(reg)
    
    B = reg(:,i); %grab column
    B = B(~isnan(B)); %remove NaNs
    
    B2 = reg(:,i+1); %grab 2nd column
    B2 = B2(~isnan(B2)); %remove NaNs
    
   [h,p1] = ttest(B(:,i),B2(:,i));
   [h,p2] = ttest(B(:,i),B2(:,i+1));
   [h,p3] = ttest(B(:,i),B2(:,i+2));
   [h,p4] = ttest(B(:,i),B2(:,i+3));
   [h,p5] = ttest(B(:,i),B2(:,i+4));
   [h,p6] = ttest(B(:,i),B2(:,i+5));
   [h,p7] = ttest(B(:,i),B2(:,i+6));

end
   
%     %pair-wise t-tests for wing regions
%     [h,p1] = ttest(regionsMedian(:,7),regionsMedian(:,16));
%     [h,p2] = ttest(regionsMedian(:,7),regionsMedian(:,10));
%     [h,p3] = ttest(regionsMedian(:,7),regionsMedian(:,13));
%     [h,p4] = ttest(regionsMedian(:,7),regionsMedian(:,1));
%     [h,p5] = ttest(regionsMedian(:,7),regionsMedian(:,19));
%     [h,p6] = ttest(regionsMedian(:,7),regionsMedian(:,14));
%     [h,p7] = ttest(regionsMedian(:,7),regionsMedian(:,24));
%     
%     [h,p8] = ttest(regionsMedian(:,16),regionsMedian(:,10));
%     [h,p9] = ttest(regionsMedian(:,16),regionsMedian(:,13));
%     [h,p10] = ttest(regionsMedian(:,16),regionsMedian(:,1));
%     [h,p11] = ttest(regionsMedian(:,16),regionsMedian(:,19));
%     [h,p12] = ttest(regionsMedian(:,16),regionsMedian(:,14));
%     [h,p13] = ttest(regionsMedian(:,16),regionsMedian(:,24));
%     
%     [h,p14] = ttest(regionsMedian(:,10),regionsMedian(:,13));
%     [h,p15] = ttest(regionsMedian(:,10),regionsMedian(:,1));
%     [h,p16] = ttest(regionsMedian(:,10),regionsMedian(:,19));
%     [h,p17] = ttest(regionsMedian(:,10),regionsMedian(:,14));
%     [h,p18] = ttest(regionsMedian(:,10),regionsMedian(:,24));
%     
%     [h,p19] = ttest(regionsMedian(:,13),regionsMedian(:,1));
%     [h,p20] = ttest(regionsMedian(:,13),regionsMedian(:,19));
%     [h,p21] = ttest(regionsMedian(:,13),regionsMedian(:,14));
%     [h,p22] = ttest(regionsMedian(:,13),regionsMedian(:,24));
%     
%     [h,p23] = ttest(regionsMedian(:,1),regionsMedian(:,19));
%     [h,p24] = ttest(regionsMedian(:,1),regionsMedian(:,14));
%     [h,p25] = ttest(regionsMedian(:,1),regionsMedian(:,24));
%     
%     [h,p26] = ttest(regionsMedian(:,19),regionsMedian(:,14));
%     [h,p27] = ttest(regionsMedian(:,19),regionsMedian(:,24));
%     
%     [h,p28] = ttest(regionsMedian(:,14),regionsMedian(:,24));
% 
%     pVals = [p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11;p12;p13;p14;p15;...
%         p16;p17;p18;p19;p20;p21;p22;p23;p24;p25;p26;p27;p28]; %store pvalues 
% 
% dlmwrite('stats_medianV_Pvals.csv',pVals)

%% BOXPLOT OF FREQUENCY WITH STATS 
close all
figure %FREQUENCIES OF HEARTS, WINGS, DH
set(gcf,'color','w') %change whole figure background to white

dat = [wingPulse(:,5), wingPulse(:,6),wingPulse(:,3),wingPulse(:,1),wingPulse(:,4),wingPulse(:,2)];
labels = {'DH','SCUT','FW Hinge','FW','HW Hinge','HW'};


%Plot
CategoricalScatterplot(dat,...
    'MedianLineWidth', 6, 'MarkerSize',50, 'WhiskerLineWidth', 4,...
    'WhiskerColor', [0.3438 0.3477 0.3477],'MedianColor',[0.4492 0.0469 0.0078],...
    'BoxColor',[0.6523 0.6602 0.6719],'BoxAlpha',1,'Labels',labels)
ylabel('Pulse frequency (Hz)','FontSize',20)
xlabel('Wing and Body Regions','FontSize',20)
ax = gca; ax.FontSize = 14; %change axes numbers to larger font


%STATS
%dat = logdat);
[~,~,stats] = anova1(dat)
figure
c = multcompare(stats)
dlmwrite('freq_wings_and_hearts_Pvals2.csv',c)

%% BOXPLOT AND STATS OF MAX VELOCITIES
%close all; 
figure
set(gcf,'color','w') %change whole figure background to white

dat = [regions(:,7),regions(:,16),regions(:,10),...
    regions(:,13),regions(:,1),...
    regions(:,4),regions(:,19),regions(:,24)];
labels = {'ABD','PRO','DH','SCUT','FW Hinge','HW Hinge','FW','HW'};

CategoricalScatterplot((dat),...
    'MedianLineWidth', 6, 'MarkerSize',50, 'WhiskerLineWidth', 4,...
    'WhiskerColor', [0.3438 0.3477 0.3477],'MedianColor',[0.4492 0.0469 0.0078],...
    'BoxColor',[0.6523 0.6602 0.6719],'BoxAlpha',1,'Labels',labels)
ylabel('Max velocity','FontSize',20)
xlabel('Wing and Body Regions','FontSize',20)
%axis([0 8 -10 80])
ax = gca; ax.FontSize = 14; %change axes numbers to larger font

%Doing stats on the log of the data -- gives it a more normal distribution


[~,~,stats] = anova1(dat);
figure
c = multcompare(stats)
dlmwrite('stats_medianV_Pvals.csv',c)

%% Plot Wing Speeds Across Span
figure
x_norm = regions(:,22)/24; x_normHW = regions(:,27)/24;
ax1 = subplot(2,1,1);
scatter(x_norm, log(regions(:,19)), 125, colorspeedF(1:end,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
colormap(colorB)
c = colorbar('Ticks', [0,1], 'TickLabels', {min(sortWingMax(:,1)), max(sortWingMax(:,1))}, 'FontSize', 14);
c.Label.String = 'Mean speed of particles (um/sec)';
c.Label.FontSize = 14;
xlabel('Normalized Wing Span')
ylabel('Log of Mean Speed (um/sec)')

ax2 = subplot(2,1,2);
scatter(x_normHW, log(regions(:,24)), 125, colorspeedH(1:end,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
colormap(colorB)
c = colorbar('Ticks', [0,1], 'TickLabels', {min(sortWingMax(:,1)), max(sortWingMax(:,1))}, 'FontSize', 14);
c.Label.String = ' Mean of particles (um/sec)';
c.Label.FontSize = 14;
xlabel('Normalized Wing Span')
ylabel('Log of Mean Speed (um/sec)')

axis([ax1 ax2],[0 1 2 11]) %NEEDS TO BE NORMALIZED FIRST
%% Plot DISTANCE TRAVELED against SPAN
% figure
% x_norm = regions(:,22)/24; 
% ax1 = subplot(2,1,1);
% scatter(x_norm, (wingDist(:,1)), 125, colorDistF(1:end,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
% colormap(colorB)
% c = colorbar('Ticks', [0,1], 'TickLabels', {min(sortWingDist(:,1)), max(sortWingDist(:,1))}, 'FontSize', 14);
% c.Label.String = 'Distance traveled by particles (um)';
% c.Label.FontSize = 14;
% xlabel('Normalized Wing Span')
% ylabel('Distance (um)')
% 
% ax2 = subplot(2,1,2);
% x_normHW = regions(:,27)/24; %XY values
% scatter(x_normHW, (wingDist(:,3)), 125, colorDistH(1:end,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
% colormap(colorB)
% c = colorbar('Ticks', [0,1], 'TickLabels', {min(sortWingDist(:,1)), max(sortWingDist(:,1))}, 'FontSize', 14);
% c.Label.String = 'Distance traveled by particles (um)';
% c.Label.FontSize = 14;
% xlabel('Normalized Wing Span')
% ylabel('Distance (um)')
% 
% axis([ax1 ax2],[0 1 -10 150]) %NEEDS TO BE NORMALIZED FIRST
%% BOX PLOTS - FOR REGIONS
figure
category = L(:,2);

logR = log(R(:,8)); %MAX speeds
h = scatterhist(X(:,1),logR,'Group',category,'Kernel','on','LineWidth',4,...
    'Marker','od+s^','MarkerSize',10);
hold on
clr = get(h(1),'colororder');
boxplot(h(2),X(:,1),category,'orientation','horizontal',...
    'label',{'','','','','','','',''},'color',clr);
boxplot(h(3),logR,category,'orientation','horizontal',...
    'label',{'','','','','','','',''},'color',clr);
set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
axis(h(1),'auto');  % Sync axes
hold off;

%% DISTANCES
% close all
% figure %total dist, start to finish
% nbins = 20;
% histogram(wingDist(:,1),nbins)
% hold on
% histogram(wingDist(:,3),nbins)
% 
% figure %cumulative total dist, start to finish
% nbins = 20;
% histogram(wingDist(:,2),nbins)
% hold on
% histogram(wingDist(:,4),nbins)
% 
% figure
% nbins = 10;
% histogram(wingDist(:,2),nbins)
% hold on
% histogram(wingDist(:,4),nbins)
% histogram(wingDist(:,1),nbins)
% hold on
% histogram(wingDist(:,3),nbins)
% 
% figure
% labels = {'FW','HW'};
% boxplot([wingDist(:,1), wingDist(:,3)],'Labels',labels,...
%     'OutlierSize',6)
% ylabel('Total distance, start to finish (um)')
% xlabel('Regions')
% 
% figure
% labels = {'FW','HW','FW-cumulative', 'HW-cumulative'};
% boxplot([wingDist(:,1), wingDist(:,3),wingDist(:,2),wingDist(:,4)],'Labels',labels,...
%     'OutlierSize',6)
% ylabel('Cumulative distance (um)')
% xlabel('Regions')
% 
% 
% figure
% labels = {'FW','HW'};
% w1 = log((wingDist(:,2)./wingDist(:,1)));
% w2 = log((wingDist(:,4)./wingDist(:,3)));
% boxplot([w1,w2],'Labels',labels,...
%     'OutlierSize',6)
% ylabel('Tot-distance as a fraction of Cumulative-Dist')
% xlabel('Regions')






%% STATS FOR WING REGIONS
close all

% %NET MOVEMENT
% figure %wing regions
% labels = {'FW','HW','LE','LE-HW','MEM','MEM-HW','TIP','TIP-HW','LATT','LATT-HW','TE','TE-HW'};
% boxplot([wingMEFW(:,1),wingMEHW(:,1),wingRegMe(:,1),wingRegMe(:,6),...
%     wingRegMe(:,2),wingRegMe(:,7),wingRegMe(:,3),wingRegMe(:,8),wingRegMe(:,4),wingRegMe(:,9),...
%     wingRegMe(:,5),wingRegMe(:,10)],'Labels',labels,...
%     'OutlierSize',6)
% ylabel('Net movement')
% xlabel('FW/HW Wing Regions')
% 
% %STATS FOR SPEED
% dat = [wingMEFW(:,1),wingMEHW(:,1),wingRegMe(:,1),wingRegMe(:,6),...
%     wingRegMe(:,2),wingRegMe(:,7),wingRegMe(:,3),wingRegMe(:,8),wingRegMe(:,4),wingRegMe(:,9),...
%     wingRegMe(:,5),wingRegMe(:,10)];
% 
% [~,~,stats] = anova1(dat);
% c = multcompare(stats);
% dlmwrite('net_movement_all_wing_regions_Pvals.csv',c)

% %PATH SINU
% figure %wing regions
% labels = {'FW','HW','LE','LE-HW','MEM','MEM-HW','TIP','TIP-HW','LATT','LATT-HW','TE','TE-HW'};
% boxplot(log([wingDist(:,1),wingDist(:,2),wingRegDIST(:,1),wingRegDIST(:,6),...
%     wingRegDIST(:,2),wingRegDIST(:,7),wingRegDIST(:,3),wingRegDIST(:,8),wingRegDIST(:,4),wingRegDIST(:,9),...
%     wingRegDIST(:,5),wingRegDIST(:,10)]),'Labels',labels,...
%     'OutlierSize',6)
% ylabel('Path sinuosity')
% xlabel('FW/HW Wing Regions')
% 
% %STATS FOR SPEED
% dat = [wingDist(:,1),wingDist(:,2),wingRegDIST(:,1),wingRegDIST(:,6),...
%     wingRegDIST(:,2),wingRegDIST(:,7),wingRegDIST(:,3),wingRegDIST(:,8),wingRegDIST(:,4),wingRegDIST(:,9),...
%     wingRegDIST(:,5),wingRegDIST(:,10)];
% 
% [~,~,stats] = anova1(dat);
% c = multcompare(stats);
% dlmwrite('path_sinu_all_wing_regions_Pvals.csv',c)

%FOR SUPPLEMENTAL FIGURE ON 19 OCTOBER 2022
%MEAN SPEED - wingReg is set to speedvar 7 right now
figure %wing regions
subplot(1,3,1)
labels = {'LE','LE-HW','MEM','MEM-HW','TIP','TIP-HW','LATT','LATT-HW','TE','TE-HW'};
% boxplot([regions(:,19),regions(:,24),wingReg(:,1),wingReg(:,6),...
%     wingReg(:,2),wingReg(:,7),wingReg(:,3),wingReg(:,8),wingReg(:,4),wingReg(:,9),...
%     wingReg(:,5),wingReg(:,10),],'Labels',labels,...
%     'OutlierSize',6)

CategoricalScatterplot([wingReg(:,1),wingReg(:,6),...
    wingReg(:,2),wingReg(:,7),wingReg(:,3),wingReg(:,8),wingReg(:,4),wingReg(:,9),...
    wingReg(:,5),wingReg(:,10)],'Labels',labels)
ylabel('Mean Velocity')
xlabel('FW/HW Wing Regions')
axis([0 11 0 1500])

subplot(1,3,2)
labels = {'LE','LE-HW','MEM','MEM-HW','TIP','TIP-HW','LATT','LATT-HW','TE','TE-HW'};
% boxplot([regions(:,19),regions(:,24),wingReg(:,1),wingReg(:,6),...
%     wingReg(:,2),wingReg(:,7),wingReg(:,3),wingReg(:,8),wingReg(:,4),wingReg(:,9),...
%     wingReg(:,5),wingReg(:,10),],'Labels',labels,...
%     'OutlierSize',6)

CategoricalScatterplot([wingRegMed(:,1),wingRegMed(:,6),...
    wingRegMed(:,2),wingRegMed(:,7),wingRegMed(:,3),wingRegMed(:,8),wingRegMed(:,4),wingRegMed(:,9),...
    wingRegMed(:,5),wingRegMed(:,10)],'Labels',labels)
ylabel('Median Velocity')
xlabel('FW/HW Wing Regions')
axis([0 11 0 1500])

subplot(1,3,3)
labels = {'LE','LE-HW','MEM','MEM-HW','TIP','TIP-HW','LATT','LATT-HW','TE','TE-HW'};
% boxplot([regions(:,19),regions(:,24),wingReg(:,1),wingReg(:,6),...
%     wingReg(:,2),wingReg(:,7),wingReg(:,3),wingReg(:,8),wingReg(:,4),wingReg(:,9),...
%     wingReg(:,5),wingReg(:,10),],'Labels',labels,...
%     'OutlierSize',6)

CategoricalScatterplot([wingRegMax(:,1),wingRegMax(:,6),...
    wingRegMax(:,2),wingRegMax(:,7),wingRegMax(:,3),wingRegMax(:,8),wingRegMax(:,4),wingRegMax(:,9),...
    wingRegMax(:,5),wingRegMax(:,10)],'Labels',labels)
ylabel('Max Velocity')
xlabel('FW/HW Wing Regions')



% %STATS FOR SPEED
% dat = [regions(:,19),regions(:,24),wingReg(:,1),wingReg(:,6),...
%     wingReg(:,2),wingReg(:,7),wingReg(:,3),wingReg(:,8),wingReg(:,4),wingReg(:,9),...
%     wingReg(:,5),wingReg(:,10),];
% 
% [~,~,stats] = anova1(dat);
% c = multcompare(stats);
% dlmwrite('max_speeds_all_wing_regions_Pvals_DEC2019.csv',c)

% WING FREQUENCIES
% figure %wing regions
% labels = {'FW','HW','LE','LE-HW','MEM','MEM-HW','TIP','TIP-HW','LATT','LATT-HW','TE','TE-HW'};
% boxplot(([wingPulse(:,1),wingPulse(:,2),wingRegFreq(:,1),wingRegFreq(:,6),...
%     wingRegFreq(:,2),wingRegFreq(:,7),wingRegFreq(:,3),wingRegFreq(:,8),wingRegFreq(:,4),wingRegFreq(:,9),...
%     wingRegFreq(:,5),wingRegFreq(:,10),]),'Labels',labels,...
%     'OutlierSize',6)
% ylabel('Frequencies')
% xlabel('FW/HW Wing Regions')
% 
% STATS FOR FREQ
% dat = [wingPulse(:,1),wingPulse(:,2),wingRegFreq(:,1),wingRegFreq(:,6),...
%     wingRegFreq(:,2),wingRegFreq(:,7),wingRegFreq(:,3),wingRegFreq(:,8),wingRegFreq(:,4),wingRegFreq(:,9),...
%     wingRegFreq(:,5),wingRegFreq(:,10),];
% 
% [~,~,stats] = anova1(dat);
% c = multcompare(stats);
% dlmwrite('freq_all_wing_regions_Pvals.csv',c)

%% Plot wing speeds based on x,y coordinates
% fig = figure;
% set(fig, 'position', [136 46 1305 733])
% %FOREWING: Plots X,Y coords of speeds -- where speed is color of the CIRCLE
% ax1 = subplot(2,1,1);
% x_norm = regions(:,22)/24; y_norm = regions(:,23)/6;
% 
% scatter(x_norm, y_norm, 200, colorMEFW(1:end,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
% colormap(colorB)
% c = colorbar('Ticks', [0,1], 'TickLabels', {'Slow', 'Fast'}, 'FontSize', 14);
% c.Label.String = 'Mean speed of particles (um/sec)';
% c.Label.FontSize = 14;
% xlabel('Normalized Wing Span')
% ylabel('Normalized Wing Chord')
% hold on
% 
% %HINDWING: Plots X,Y coords of speeds -- where speed is color of the CIRCLE
% ax2 = subplot(2,1,2);
% x_norm = regions(:,27)/24; y_norm = regions(:,28)/13;
% 
% scatter(x_norm, y_norm, 200, colorMEHW(1:end,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
% colormap(colorB)
% c = colorbar('Ticks', [0,1], 'TickLabels', {'Slow', 'Fast'}, 'FontSize', 14);
% c.Label.String = 'Mean speed of particles (um/sec)';
% c.Label.FontSize = 14;
% xlabel('Normalized Wing Span')
% ylabel('Normalized Wing Chord')
% %title('Range of speeds within S. americana hindwing')
% 
% axis([ax1 ax2],[0 1 0 1])
%% DEC 2019 
% Plot TOTAL DISTANCE/SPEED/FREQ on x,y coordinates
% EDITING TO HAVE A NEW COLOR PALETTE, ORIGINAL BELOW
% 
% fig1 = figure;
% set(fig1, 'position', [1 552 1440 223])
% FOREWING: Plots X,Y coords of dist -- where dist is color of the CIRCLE
% ax1 = subplot(2,1,1);
% normalization is based on how the wings were initially divided 24
% divisions along the span for FW and HW. 6 divisions on the FW-y-axis and
% 13 divisions on the HW-y-axis
% x_norm = fwEDIT_outX/24; y_norm = fwEDIT_outY/6; 
% 
% 
% minFW = min(regions(:,19)); minHW = min(regions(:,24));
% maxFW = max(regions(:,19)); maxHW = max(regions(:,24));
% 
% minF = min([minFW, minHW]); %get min for both FW/HW
% maxF = max([maxFW, maxHW]); %get max for both FW/HW
% 
% scatter(x_norm, y_norm, 250, colorspeedF_EDIT(1:end,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
% colormap(flipcolorD)
% c = colorbar('Ticks', [0,1], 'TickLabels', {minF,maxF}, 'FontSize', 14);
% c.Label.String = 'Velocity (um/s)';
% c.Label.FontSize = 18;
% xlabel('Normalized Wing Span','FontSize',18)
% ylabel('Normalized Wing Chord','FontSize',18)
% axis([-0.2 1.2 -0.2 1.2])
% hold on

% %% DEC EDITS 2019
% %HINDWING: Plots X,Y coords of DISTANCE -- where DISTANCE is color of the CIRCLE
% %ax2 = subplot(2,1,2);
% fig2 = figure;
% set(fig2, 'position', [136 46 1305 733])
% %normalization is based on how the wings were initially divided 24
% %divisions along the span for FW and HW. 6 divisions on the FW-y-axis and
% %13 divisions on the HW-y-axis
% x_norm = hwEDIT_outX/24; y_norm = hwEDIT_outY/13;
% 
% scatter(x_norm, y_norm, 250, colorspeedH_EDIT(1:end,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
% colormap(flipcolorD)
% c = colorbar('Ticks', [0,1], 'TickLabels', {minF,maxF}, 'FontSize', 14);
% c.Label.String = 'Velocity (um/s) ';
% c.Label.FontSize = 18;
% xlabel('Normalized Wing Span','FontSize',18)
% ylabel('Normalized Wing Chord','FontSize',18)
% axis([-0.2 1.2 -0.2 1.2])
% 
% %axis([ax1 ax2],[0 1 0 1])

%% ORIGINAL: Plot TOTAL DISTANCE/SPEED/FREQ on x,y coordinates


fig1 = figure;
set(fig1, 'position', [1 552 1440 223])
%FOREWING: Plots X,Y coords of dist -- where dist is color of the CIRCLE
%ax1 = subplot(2,1,1);
%normalization is based on how the wings were initially divided 24
%divisions along the span for FW and HW. 6 divisions on the FW-y-axis and
%13 divisions on the HW-y-axis
x_norm = regions(:,22)/24; y_norm = regions(:,23)/6; 


minFW = min(regions(:,19)); minHW = min(regions(:,24));
maxFW = max(regions(:,19)); maxHW = max(regions(:,24));

minF = min([minFW, minHW]); %get min for both FW/HW
maxF = max([maxFW, maxHW]); %get max for both FW/HW

scatter(x_norm, y_norm, 200, colorFreqPFW(1:end,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
colormap(colorF)
c = colorbar('Ticks', [0,1], 'TickLabels', {minF,maxF}, 'FontSize', 14);
c.Label.String = 'Velocity (um/s)';
c.Label.FontSize = 18;
xlabel('Normalized Wing Span','FontSize',18)
ylabel('Normalized Wing Chord','FontSize',18)
axis([-0.2 1.2 -0.2 1.2])
hold on

%HINDWING: Plots X,Y coords of DISTANCE -- where DISTANCE is color of the CIRCLE
%ax2 = subplot(2,1,2);
fig2 = figure;
set(fig2, 'position', [136 46 1305 733])
%normalization is based on how the wings were initially divided 24
%divisions along the span for FW and HW. 6 divisions on the FW-y-axis and
%13 divisions on the HW-y-axis
x_norm = regions(:,27)/24; y_norm = regions(:,28)/13;

scatter(x_norm, y_norm, 200, colorFreqPHW(1:end,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
colormap(colorF)
c = colorbar('Ticks', [0,1], 'TickLabels', {minF,maxF}, 'FontSize', 14);
c.Label.String = 'Velocity (um/s) ';
c.Label.FontSize = 18;
xlabel('Normalized Wing Span','FontSize',18)
ylabel('Normalized Wing Chord','FontSize',18)
axis([-0.2 1.2 -0.2 1.2])

%axis([ax1 ax2],[0 1 0 1])
%% Quiver plots -- MEAN VELOCITY
%arrow color
c = [128/255, 130/255, 133/255]; %gray

k = find(R(:,6) == 0); %index for L wings
k2 = find(R(:,6) == 1); %index for R wings

%scalefactors and normalization
xmax = 1.2; xmin = -0.2; ymax = 1.2; ymin = -0.2;
x_span = xmax - xmin; y_span = ymax - ymin;
d = 0.01;
s1_fw = 6*(d./sqrt((regions(:,20)/x_span).^2+(regions(:,21)/y_span).^2));
s2_hw = 6*(d./sqrt((regions(:,25)/x_span).^2+(regions(:,26)/y_span).^2));

%left, fw wings --  based on index k = 0 (L)
x_L = regions(k,22)/24; %grab positions of left wings
y_L = regions(k,23)/6; %flip the y axis for left wings

speedVarL = regions(k,19); 
u_L = -s1_fw(k).*regions(k,20); %x and y speeds
v_L = s1_fw(k).*regions(k,21); 

%left, hw wings  --  based on index k = 0 (L)
x_L_HW = regions(k,27)/24; %grab posionts of left wings
y_L_HW = regions(k,28)/13; %flip the y axis for left wings

speedVarLHW = regions(k,24);
uL_HW = -s2_hw(k).*regions(k,25); %x and y speeds
vL_HW = s2_hw(k).*regions(k,26); 

%right fw wings  --  based on index k2 = 1 (R)
x_R =  regions(k2,22)/24; %grab positions of right wings
y_R =  regions(k2,23)/6;

speedVarR = regions(k2,19); %velocity FW values (RIGHT)
u_R = s1_fw(k2).*regions(k2,20); %x and y vectors, (u,v)
v_R = s1_fw(k2).*regions(k2,21);

%right hw wings  --  based on index k2 = 1 (R)
x_R_HW =  regions(k2,27)/24; %grab positions of right wings
y_R_HW =  regions(k2,28)/13; 

speedVarRHW = regions(k2,24); %velocity HW values (RIGHT)
uR_HW = s2_hw(k2).*regions(k2,25); %x and y speeds
vR_HW = s2_hw(k2).*regions(k2,26);

%% QUIVER EDITS DECEMBER 2019
%arrow color
c = [128/255, 130/255, 133/255]; %gray

k = find(R(:,6) == 0); %index for L wings
k2 = find(R(:,6) == 1); %index for R wings

%scalefactors and normalization
xmax = 1.2; xmin = -0.2; ymax = 1.2; ymin = -0.2;
x_span = xmax - xmin; y_span = ymax - ymin;
d = 0.01;
s1_fw = 6*(d./sqrt((regions(:,20)/x_span).^2+(regions(:,21)/y_span).^2));
%s1_fw = s1_fw(~isnan(s1_fw));

s2_hw = 6*(d./sqrt((regions(:,25)/x_span).^2+(regions(:,26)/y_span).^2));
%s2_hw = s2_hw(~isnan(s2_hw));

%left, fw wings --  based on index k = 0 (L)
%hwEDIT_outX = hwEDIT_X(~isnan(hwEDIT_X));
x_L = regions(k,22)/24; %grab positions of left wings
x_L = x_L(~isnan(x_L)); %186

y_L = regions(k,23)/6; %flip the y axis for left wings
y_L = y_L(~isnan(y_L)); %186

speedVarL = regions(k,19); 
speedVarL = speedVarL(~isnan(speedVarL)); %186

u_L = -s1_fw(k).*regions(k,20); %x and y speeds
u_L = u_L(~isnan(u_L));

v_L = s1_fw(k).*regions(k,21); 
v_L = v_L(~isnan(v_L));

%left, hw wings  --  based on index k = 0 (L)
x_L_HW = regions(k,27)/24; %grab posionts of left wings
x_L_HW = x_L_HW(~isnan(x_L_HW));
y_L_HW = regions(k,28)/13; %flip the y axis for left wings
y_L_HW = y_L_HW(~isnan(y_L_HW));

speedVarLHW = regions(k,24);
speedVarLHW = speedVarLHW(~isnan(speedVarLHW));
uL_HW = -s2_hw(k).*regions(k,25); %x and y speeds
uL_HW = uL_HW(~isnan(uL_HW));
vL_HW = s2_hw(k).*regions(k,26); 
vL_HW = vL_HW(~isnan(vL_HW));

%right fw wings  --  based on index k2 = 1 (R)
x_R =  regions(k2,22)/24; %grab positions of right wings
x_R = x_R(~isnan(x_R));
y_R =  regions(k2,23)/6;
y_R = y_R(~isnan(y_R));

speedVarR = regions(k2,19); %velocity FW values (RIGHT)
speedVarR = speedVarR(~isnan(speedVarR));
u_R = s1_fw(k2).*regions(k2,20); %x and y vectors, (u,v)
u_R = u_R(~isnan(u_R));
v_R = s1_fw(k2).*regions(k2,21);
v_R = v_R(~isnan(v_R));

%right hw wings  --  based on index k2 = 1 (R)
x_R_HW =  regions(k2,27)/24; %grab positions of right wings
x_R_HW = x_R_HW(~isnan(x_R_HW));
y_R_HW =  regions(k2,28)/13; 
y_R_HW = y_R_HW(~isnan(y_R_HW));

speedVarRHW = regions(k2,24); %velocity HW values (RIGHT)
speedVarRHW = speedVarRHW(~isnan(speedVarRHW));
uR_HW = s2_hw(k2).*regions(k2,25); %x and y speeds
uR_HW = uR_HW(~isnan(uR_HW));
vR_HW = s2_hw(k2).*regions(k2,26);
vR_HW = vR_HW(~isnan(vR_HW));


%%
%Make quiver plot of speeds on the wings
fig1 = figure;
set(fig1, 'position', [1 552 1440 223])

%QUIV PROPS
line_w = 2; n = 'off'; auto = 4; maxHS = 2; ah = 'on';

%GOAL: plot RIGHT FOREWINGS
for i = 1:length(x_R)
    
    indx = find(speedVarR(i) <= sortcombo(:,1));
    if ~isnan(indx)
        if isempty(indx)
            indx = find(speedVarR(i)<=sortcombo(:,1));
        end
        indx = round(indx(1)/10); %since high #s, divide and round up
        
        if indx == 0 
            indx = 1;
        end
        quiver(x_R(i), y_R(i), u_R(i), v_R(i),'Color',flipcolorB(indx,:),'LineWidth',line_w,...
            'ShowArrowHead',ah,'MaxHeadSize', maxHS,'AutoScaleFactor',auto, 'AutoScale',n);
        
    end
    hold on
end
%LEFT FOREWINGS
for i = 1:length(x_L)
    indx = find(speedVarL(i) <= sortcombo(:,1));
    if ~isnan(indx)
        if isempty(indx)
            indx = find(speedVarL(i)<=sortcombo(:,1));
        end
        indx = round(indx(1)/10); %since high #s, divide and round up
        
        if indx == 0 
            indx = 1;
        end
        quiver(x_L(i), y_L(i), u_L(i), v_L(i),'Color',flipcolorB(indx,:),'LineWidth',line_w,...
            'ShowArrowHead',ah,'MaxHeadSize', maxHS,'AutoScaleFactor',auto, 'AutoScale',n);
        
    end
    hold on
end
fig1=figure(1);
fig1.Renderer='Painters';

%ADD VELOCITY DOTS - forewing
%x_norm = regions(:,22)/24; y_norm = regions(:,23)/6;
x_norm = fwEDIT_outX/24; y_norm = fwEDIT_outY/6; %follow the code, this is calculated above

scatter(x_norm, y_norm, 75, colorspeedF_EDIT(1:end,:), 'filled')
colormap(flipcolorB)
f = colorbar('Ticks', [0,1], 'TickLabels', {minF, maxF}, 'FontSize', 14); 
f.Label.String = 'Min speed of particles (um/sec)'; f.Label.FontSize = 14;
xlabel('Normalized Wing Span')
ylabel('Normalized Wing Chord')
axis([-0.2 1.2 -0.2 1.2])
hold off

fig2 = figure;
set(fig2, 'position', [136 46 1305 733])
%subplot(2,1,2)
%plot R and L HWs
for i = 1:length(x_R_HW)
    indx = find(speedVarRHW(i) <= sortcombo(:,1));
    if ~isnan(indx)
        if isempty(indx)
            indx = find(speedVarRHW(i)<=sortcombo(:,1));
        end
        indx = round(indx(1)/10); %since high #s, divide and round up
 
        if indx == 0 
            indx = 1;
        end
        quiver(x_R_HW(i), y_R_HW(i), uR_HW(i), vR_HW(i),'Color',flipcolorB(indx,:),'LineWidth',line_w,...
            'ShowArrowHead',ah,'MaxHeadSize', maxHS,'AutoScaleFactor',auto, 'AutoScale',n);
    end
    hold on 
end 
for i = 1:length(x_L_HW)
    indx = find(speedVarLHW(i) <= sortcombo(:,1));
    if ~isnan(indx)
        if isempty(indx)
            indx = find(speedVarLHW(i)<=sortcombo(:,1));
        end
        indx = round(indx(1)/10); %since high #s, divide and round u
        
        if indx == 0 
            indx = 1;
        end
            quiver(x_L_HW(i), y_L_HW(i), uL_HW(i), vL_HW(i),'Color',flipcolorB(indx,:),'LineWidth',line_w,...
                'ShowArrowHead',ah,'MaxHeadSize', maxHS,'AutoScaleFactor',auto, 'AutoScale',n);
    end
    hold on
end

%ADD VELOCITY DOTS - hindwing 
%x_norm = regions(:,27)/24; y_norm = regions(:,28)/13;
x_norm = hwEDIT_outX/24; y_norm = hwEDIT_outY/13; % follow the code - calculated modules above

scatter(x_norm, y_norm, 75, colorspeedH_EDIT(1:end,:), 'filled') %colorspeedH_EDIT calc'd above
colormap(flipcolorB) %colorbar of 34 vals
f = colorbar('Ticks', [0,1], 'TickLabels', {minF, maxF}, 'FontSize', 14);
f.Label.String = 'Min speed of particles (um/sec)'; f.Label.FontSize = 14;
xlabel('Normalized Wing Span')
ylabel('Normalized Wing Chord')
axis([-0.2 1.2 -0.2 1.2])
fig2=figure(2);
fig2.Renderer='Painters';

hold off
%%
% %%
% %plot R and L FWs
% line_w = 2; n = 'on'; auto = 1; maxHS = 2;
% 
% %GOAL: plot FOREWINGS
% for i = 1:length(x_R)
%     nearest_index = @(vector,element) find(abs(element-vector) == min(abs(element-vector)));
%     vector = sortMax(:,1);
%     element = meanR(i);
%     indx_F = nearest_index(vector,element);
%     %value = vector(index);
% 
%     %f = find(abs(meanR(i)) == sortMax(:,1)); %find indx of speeds
%     indx = round(indx_F/10); %since high #s, divide and round up
%         
%     if ~isempty(indx)
%         if indx(1) ==0
%             indx = 1;
%         end
%         quiver(x_R(i), y_R(i), u_R(i), v_R(i),'Color',colorB2(indx(1),:),'LineWidth',line_w,...
%             'MaxHeadSize', maxHS,'AutoScaleFactor',auto, 'AutoScale',n);
%     end
%     hold on
% end
% %LEFT FOREWINGS
% for i = 1:length(x_L)
%     nearest_index = @(vector,element) find(abs(element-vector) == min(abs(element-vector)));
%     vector = sortMax(:,1);
%     element = meanL(i);
%     indx_F = nearest_index(vector,element);
%     %value = vector(index);
% 
%     %f = find(abs(meanR(i)) == sortMax(:,1)); %find indx of speeds
%     indx = round(indx_F/10); %since high #s, divide and round up
%     
%     
%     if ~isempty(indx)
%         if indx(1)==0
%             indx = 1;
%         end
%         quiver(x_L(i), y_L(i), u_L(i), v_L(i),'Color',colorB2(indx(1),:),'LineWidth',line_w,...
%             'MaxHeadSize', maxHS,'AutoScaleFactor',auto, 'AutoScale',n);
%     end
%     hold on
% end
% fig1=figure(1);
% fig1.Renderer='Painters';
% 
% %add velocity dots
% x_norm = regions(:,22)/24; y_norm = regions(:,23)/6;
% scatter(x_norm, y_norm, 75, colorspeedmeanF(1:end,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
% colormap(colorB)
% f = colorbar('Ticks', [0,1], 'TickLabels', {min(sortMax(:,1)), max(sortMax(:,1))}, 'FontSize', 14); 
% f.Label.String = 'Min speed of particles (um/sec)'; f.Label.FontSize = 14;
% xlabel('Normalized Wing Span')
% ylabel('Normalized Wing Chord')
% axis([-0.2 1.2 -0.2 1.2])
% hold off
% 
% fig2 = figure;
% set(fig2, 'position', [136 46 1305 733])
% %subplot(2,1,2)
% %plot R and L HWs
% for i = 1:length(x_R_HW)
%     nearest_index = @(vector,element) find(abs(element-vector) == min(abs(element-vector)));
%     vector = sortMax(:,1);
%     element = meanRHW(i);
%     indx_F = nearest_index(vector,element);
%     %value = vector(index);
% 
%     %f = find(abs(meanR(i)) == sortMax(:,1)); %find indx of speeds
%     indx = round(indx_F/10); %since high #s, divide and round up
%     
%     
%     if ~isempty(indx)
%         if indx(1)==0
%             indx = 1;
%         end
%         quiver(x_R_HW(i), y_R_HW(i), uR_HW(i), vR_HW(i),'Color',colorB2(indx(1),:),'LineWidth',line_w,...
%             'MaxHeadSize', maxHS,'AutoScaleFactor',auto, 'AutoScale',n);
%     end
%     hold on 
% end 
% for i = 1:length(x_L_HW)
%     nearest_index = @(vector,element) find(abs(element-vector) == min(abs(element-vector)));
%     vector = sortMax(:,1);
%     element = meanLHW(i);
%     indx_F = nearest_index(vector,element);
%     %value = vector(index);
% 
%     %f = find(abs(meanR(i)) == sortMax(:,1)); %find indx of speeds
%     indx = round(indx_F/10); %since high #s, divide and round up
%     
%     
%     if ~isempty(indx)
%         if indx(1)==0
%             indx = 1;
%         end
%         quiver(x_L_HW(i), y_L_HW(i), uL_HW(i), vL_HW(i),'Color',colorB2(indx(1),:),'LineWidth',line_w,...
%             'MaxHeadSize', maxHS,'AutoScaleFactor',auto, 'AutoScale',n);
%     end
%     hold on
% end
% 
% %add velocity dots
% x_norm = regions(:,27)/24; y_norm = regions(:,28)/13;
% scatter(x_norm, y_norm, 75, colorspeedmeanH(1:end,:), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5])
% colormap(colorB)
% f = colorbar('Ticks', [0,1], 'TickLabels', {min(sortMax(:,1)), max(sortMax(:,1))}, 'FontSize', 14);
% f.Label.String = 'Min speed of particles (um/sec)'; f.Label.FontSize = 14;
% xlabel('Normalized Wing Span')
% ylabel('Normalized Wing Chord')
% axis([-0.2 1.2 -0.2 1.2])
% fig2=figure(2);
% fig2.Renderer='Painters';
% 
% hold off




