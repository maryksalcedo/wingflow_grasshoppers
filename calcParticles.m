%made to return quiver values
function [particleMean,numDig,particleMax,partMed,numPks,cum_dist,alpha_1,re_1] = calcParticles(xyPTS,voxSize,timeInt)
    A = dlmread(xyPTS, ',', 1,0); %this should read in all of the xypts.csv
    numPts = round(size(A,2)); %round to an even num of digitized points (aka the columns)

    %OLD : A_diff = (diff(A)*voxSize)/timeInt; %take diff, change pixels to um
    A_diff = (diff(A)*voxSize); %NEW: take diff, change PIXELS ARE NOW MICRONS
    
    %set up empties
    dist_path = zeros(length(A_diff),numPts); %predefine total distance
    vMag = zeros(length(A_diff),numPts); %predefine velocity matrix
    x_diff = zeros(length(A_diff),numPts/2); %predefine x-diff
    y_diff = zeros(length(A_diff),numPts/2); %predefine y-diff
    
    

    k = 1;
    for i = 1:2:numPts %Loop through all of columns of x,y pairs in the data xypts matrix
        x_diff(:,i) = A_diff(:,i); % column of x-value
        y_diff(:,i) = A_diff(:,i+1); %column of y-values
        
        %x_pos(:,i) = A_pos(:,i); %x-value
        %y_pos(:,i) = A_pos(:,i+1); %y-values
        
        %calculate distance traveled
        %this is incorrect - it does not subtract
        %OLD dist_path = sqrt(x_pos(:,i).^2 + y_pos(:,i).^2);
      
        dist_path = sqrt(x_diff(:,i).^2 + y_diff(:,i).^2); %calculates distance
        cum_dist(k,1) = sum(dist_path,'omitnan'); %cumulative distance traveled
        
%         list = ~isnan(dist_path); %get rid of NaNs
%         values_cum = dist_path(list,i); 
%         indx = find(~isnan(x_pos(:,i))); %find 1st index when not an NaN
%         x1 = x_pos(indx(1),i); x2 = x_pos(indx(end),i); y1 = y_pos(indx(1),i); y2 = y_pos(indx(end),i);    
%         %Net_dist does not take into account vein structure (necessarily)
%         net_dist(k,1) = sqrt(diff([x1;x2]).^2 + diff([y1;y2]).^2); %overall, start-to-end distance
        
       net_dist = 0;
        %VELOCITY: instantaneous velocity -- distance calc with time
        b = x_diff(:,i); numDig(k,1) = length(b(~isnan(b))); %number of digitized frames
        time = numDig(k,1)*timeInt; %Time is in seconds
        %vMag(:,i) = sqrt(x_diff(:,i).^2 + y_diff(:,i).^2)/time; %microns/sec
           
        vMag(:,i) = sqrt(x_diff(:,i).^2 + y_diff(:,i).^2)/timeInt; %microns/sec
        
        %FIND PEAKS:
        xdat = abs(movmean(x_diff(:,i),10,'omitnan')); %over-smooth and absolute value
        norm_xdat = xdat/max(xdat); %normalize by dividing by max
        [pks,~] = findpeaks(norm_xdat,'MinPeakProminence', 0.3); %0.3 picks up most peaks adjust as needed
        numPks(k,1) = length(pks);
        

        %Womersley = Alpha = L*(omega*density/dynamic viscosity)^(1/2)
        %This Womersley will leave out radius of pipe until the code sorts
        %per region
        omega = numPks(k,1)/time; %frequency (times per second)
        dyn_vis = 0.0010518; %Pa*sec of water at 18 C / 64 F
        density = 997; %kg/m^3 of water
        alpha_1(k,1) = ((omega*density)/(dyn_vis)).^(0.5); %part 1 of Womersley
        
        
        particleMean(k,1) = mean(movmean(vMag(:,i),2,'omitnan'),'omitnan');
        particleMax(k,1) = max(movmean(vMag(:,i),2,'omitnan'));
        partMed(k,1) = median(vMag(:,i),'omitnan');
        
        %reynolds number = rho*V*radius/dyn viscosity
        %this will be Re part 1
        re_1(k,1) = (particleMax(k,1)*(1*10^(-6))*density)/dyn_vis; %convert microns to meters (10^-6)
    
        k = k+1;
    end
end