%This sript is made by ROB18_Gr465 AAU 

%This sript performs the multilateration method on the recordings, recorded
%from play_rec.

% This sript uses two add-Ons from the matlab add-on library. The first
% add-on is made by Douglas Schwarz and is called 
% Fast and Robust Curve Intersections. It is uses to find all intersection
% between two function.
% The second add-on is made by Yann Marcon and is called Distance-based 
% clustering of a set of XY coordinates. This function finds clusters in a
% set of spatial points expressed in
% XY coordinates. The clustering is based on the distance between the
% points and it does not require the number of clusters to be known
% beforehand. Each point is clustered with the closest neighbouring point
% if the distance between the two points is shorter than the user-defined
% threshold.

%links for both add-Ons
%https://se.mathworks.com/matlabcentral/fileexchange/56150-distance-based-
%clustering-of-a-set-of-xy-coordinates

%https://se.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-
%curve-intersections

function [] = New_test(rec)
close all
%% Setup %%

v_sound = 343.2;        %Speed of sound
fs = 96000;             %Sampling freqency

%Splitting the data up into 3 for each microphones
mic_1 = -(rec(1,:));    
mic_2 = -(rec(2,:));    
mic_3 = -(rec(3,:));    
    
%% %Plot Raw data
figure(1)
subplot(3,1,1)
plot(mic_1)
title('Raw Data 1','FontSize',20)
subplot(3,1,2)
plot(mic_2)
title('Raw Data 2','FontSize',20)
subplot(3,1,3)
plot(mic_3)
title('Raw Data 3','FontSize',20)

%Making ButterWorth filter 
[d1,d2] = butter(5, 3000/(fs/2), 'high');

%Using the filter on all mic data
x1 = filtfilt(d1,d2,double(mic_1));
x2 = filtfilt(d1,d2,double(mic_2));
x3 = filtfilt(d1,d2,double(mic_3));

%% %Plot filtered data
figure(3)
subplot(3,1,1)
plot(x1)
title('Filtered Data 1','FontSize',20)
subplot(3,1,2)
plot(x2)
title('Filtered Data 2','FontSize',20)
subplot(3,1,3)
plot(x3)
title('Filtered Data 3','FontSize',20)

%% setup for parabola
%Mic position on the x-axis
mic_1x = 0.15;
mic_2x = 0;
mic_3x = -0.15;

%% Processing the data 

%Finding the maximum amplitude for all mics 
[M1,~] = max(mic_1);
[M2,~] = max(mic_2);
[M3,~] = max(mic_3);

%Find all peaks above the theashold (Max_amplitude/2)
[~, index1] = findpeaks(mic_1, 'MinPeakDistance', 1, 'MinPeakHeight', M1/2);
[~, index2] = findpeaks(mic_2, 'MinPeakDistance', 1, 'MinPeakHeight', M2/2);
[~, index3] = findpeaks(mic_3, 'MinPeakDistance', 1, 'MinPeakHeight', M3/2);

%Removing of the intial pulse
temp1 = x1;
for i = index1(1)-20:index1(1)+100
    temp1(i) = 0;
end 

%Removing of the intial pulse
temp2 = x2;
for i = index2(1)-20:index2(1)+100
    temp2(i) = 0;
end 

%Removing of the intial pulse
temp3 = x3;
for i = index3(1)-20:index3(1)+100
    temp3(i) = 0;
end 


%% %Plot filtered data without the intial pulse
figure(5)
subplot(3,1,1)
plot(temp1)
title('Removing first spike 1','FontSize',20)
subplot(3,1,2)
plot(temp2)
title('Removing first spike 2','FontSize',20)
subplot(3,1,3)
plot(temp3)
title('Removing first spike 3','FontSize',20)



%% Offset x-axis depending of mic position
x = -3:0.001:3;

x12 = -3+mic_1x/2:0.001:3+mic_1x/2;
x23 = -3+mic_3x/2:0.001:3+mic_3x/2;
x13 = -3:0.001:3;

%Calculating the distance between the microphones
D_12 = sqrt((mic_1x-mic_2x)^2);
D_23 = sqrt((mic_2x-mic_3x)^2);
D_13 = sqrt((mic_1x-mic_3x)^2);

%% Cross Correlation

%Calculating the cross correlation for all pairs of microphones
[acor_12 lag_12] = xcorr(temp1,temp2);
[acor_23 lag_23] = xcorr(temp2,temp3);
[acor_13 lag_13] = xcorr(temp1,temp3);

%Deleting everything beyond the field of view(0 to 180 degress) infront of
%the mics in lags
FOV = ceil((D_13/v_sound)*fs)
for i = 1:ceil((length(acor_12)/2))-FOV
    acor_12(i) = 0;
    acor_23(i) = 0;
    acor_13(i) = 0;
    
    acor_12(length(acor_12)-i) = 0;
    acor_23(length(acor_12)-i) = 0;
    acor_13(length(acor_12)-i) = 0;
end

%% Selecting of TDOAs

%Finding the maxium amplitude of the cross correlation
M12 = max(acor_12);
M23 = max(acor_23);
M13 = max(acor_13);

%Find all peaks above the theashold (Max_amplitude/threshold)
threshold = 2;
[~, I1] = findpeaks(acor_12, 'MinPeakDistance', 1, 'MinPeakHeight', M12/threshold);
[~, I2] = findpeaks(acor_23, 'MinPeakDistance', 1, 'MinPeakHeight', M23/threshold);
[~, I3] = findpeaks(acor_13, 'MinPeakDistance', 1, 'MinPeakHeight', M13/threshold);


%% %Plot the cross correlation with the peaks above the threshold
figure(2)
clf
findpeaks(acor_12, 'MinPeakDistance', 1, 'MinPeakHeight', M12/threshold)
hold on
findpeaks(acor_23, 'MinPeakDistance', 1, 'MinPeakHeight', M23/threshold)
hold on
findpeaks(acor_13, 'MinPeakDistance', 1, 'MinPeakHeight', M13/threshold)
xlim([8100 8280])

%% Calculating Hyperbolas

%Calculating the TDOAs for all pairs of microphones
TDOA_12 = lag_12(I1)/fs;
TDOA_23 = lag_23(I2)/fs;
TDOA_13 = lag_13(I3)/fs;

%Calculating delta_d based on the TDOAs
d_12 = v_sound*TDOA_12;
d_23 = v_sound*TDOA_23;
d_13 = v_sound*TDOA_13;

%Removing all zeros from the delta_d
d_12(d_12 == 0) = []
d_23(d_23 == 0) = []
d_13(d_13 == 0) = []

%Initialize hyperbola arrays for all mic pairs
Y12 = [];       
Y23 = [];       
Y13 = [];

%Calculation the hyperbolas based on the nummber of TDOAs
for i = 1:length(d_12)
Y12 = [Y12; 1/2*sqrt(((4*x12.^2-d_12(i)^2)*(D_12^2-d_12(i)^2))/d_12(i)^2)];
end

for i = 1:length(d_23)
Y23 = [Y23; 1/2*sqrt(((4*x23.^2-d_23(i)^2)*(D_23^2-d_23(i)^2))/d_23(i)^2)];
end

for i = 1:length(d_13)
Y13 = [Y13; 1/2*sqrt(((4*x13.^2-d_13(i)^2)*(D_13^2-d_13(i)^2))/d_13(i)^2)];
end

%% %Plot all hyperbolas
figure(4)
clf
if ~isempty(Y12)
plot(x,real(Y12), 'r');
hold on
end
if ~isempty(Y23)
plot(x,real(Y23), 'g'); 
hold on
end
if ~isempty(Y13)
plot(x,real(Y13), 'b');
hold on
end
title({'Mic setup = -0.15, 0, 0.15'; 'Source: 0.5 m, -45 degrees'})
xlim([-2 2]);ylim([-2 4]);

%% %Find intersections

%Initialize arrays for storing X and Y coordinates for all intersections
XOUT = [];
YOUT = [];

%Find all intersection between hypperbolas 12 and 23 and store the
%coordinates
if ~isempty(Y12) && ~isempty(Y23)
    XOUT12= [];
    YOUT12= [];
    [m12, ~] = size(Y12);
    [m23, ~] = size(Y23);
    for i = 1:m12
        temp1 = Y12(i,:);
        for j = 1:m23
            temp2 = Y23(j,:);
            %Using a intersection add-on called Fast and Robust Curve
            %Intersections. This finds all intersections between to
            %fuctions
            [XOUT12a, YOUT12a] = intersections(x,temp1,x,temp2,true);
            XOUT12 = [XOUT12; XOUT12a];
            YOUT12 = [YOUT12; YOUT12a];
        end
    end
    
    %Storing X and Y coordinate in array
    [~, n12]=size(XOUT12);
    for i = 1:n12
        XOUT = [XOUT; real(XOUT12(:,i))];
    end
    [~, n12]=size(YOUT12);
    for i = 1:n12
        YOUT = [YOUT; real(YOUT12(:,i))];
    end
end

%Find all intersection between hypperbolas 23 and 13 and store the
%coordinates
if ~isempty(Y23) && ~isempty(Y13)
    XOUT23= [];
    YOUT23= [];
    [m23, ~] = size(Y23);
    [m13, ~] = size(Y13);
    for i = 1:m23
            temp1 = Y23(i,:);
        for j = 1:m13
            temp2 = Y13(j,:);
            %Using a intersection add-on called Fast and Robust Curve
            %Intersections. This finds all intersections between to
            %fuctions
            [XOUT23a, YOUT23a] = intersections(x,temp1,x,temp2,true);
            XOUT23 = [XOUT23; XOUT23a];
            YOUT23 = [YOUT23; YOUT23a];
        end
    end
    
    %Storing X and Y coordinate in array
    [~, n23]=size(XOUT23);
    for i = 1:n23
        XOUT = [XOUT; real(XOUT23(:,i))];
    end
    [~, n23]=size(YOUT23);
    for i = 1:n23
        YOUT = [YOUT; real(YOUT23(:,i))];
    end
end

%Find all intersection between hypperbolas 13 and 12 and store the
%coordinates
if ~isempty(Y13) && ~isempty(Y12)
    XOUT13= [];
    YOUT13= [];
    [m13, ~] = size(Y13);
    [m12, ~] = size(Y12);
    
    for i = 1:m13
            temp1 = Y13(i,:);
        for j = 1:m12
            temp2 = Y12(j,:);
            %Using a intersection add-on called Fast and Robust Curve
            %Intersections. This finds all intersections between to
            %fuctions
            [XOUT13a, YOUT13a] = intersections(x,temp1,x,temp2,true);
            XOUT13 = [XOUT13; XOUT13a];
            YOUT13 = [YOUT13; YOUT13a];
        end
    end
    
    %Storing X and Y coordinate in array
    [~, n13]=size(XOUT13);
    for i = 1:13
        XOUT = [XOUT; real(XOUT13(:,i))];
    end
    [~, n13]=size(YOUT13);
    for i = 1:n13
        YOUT = [YOUT; real(YOUT13(:,i))];
    end
end

Coord = [];
CoordXY= [XOUT YOUT];

%Saving all intersections that have a lenght larger than 0.15m and less than 5m to
%the origin
[A,~] = size(CoordXY);
for i = 1:(A)
   if sqrt(CoordXY(i,1)^2+CoordXY(i,2)^2) > 0.15 && sqrt(CoordXY(i,1)^2+CoordXY(i,2)^2) < 5
      Coord(i,:) = CoordXY(i,:);
   end
end

%Removing all zeros from the intersections
[row,~] = find(Coord(:,1) ~= 0)
new_coords = Coord(row,:)

%% %Plot all valid intersections
scatter(new_coords(:,1), new_coords(:,2), 'filled');
title('Hyperbolas, shown positively for positive and negative lag')
xlabel('Distance [m]')
ylabel('Distance [m]')
hold on

%Plot microphones position
pos_mic = [0, 0;0.15, 0;-0.15,0]
scatter(pos_mic(:,1),pos_mic(:,2),100,'X')

%Quick fix for clustering
new_coords = [new_coords ; 0 , 0];

%% %Clustering to find objects
%Using a clustering add-on called Distance-based clustering of a set of XY
%coordinates. This is used to find 3 points that have a distance between
%them less than 0.0005m. and finds the mean for the 3 points.
[clustersCentroids, clustersGeoMedians, clustersXY] = clusterXYpoints(new_coords, 0.0005, 3, 'point', 'merge');
clustersCentroids

%% %Plot the clusterCentroids
if ~isempty(clustersCentroids)
    abc = scatter(clustersCentroids(:,1),clustersCentroids(:,2), 200, 'filled', 'R')
    abc.MarkerEdgeColor = 'k';
    ylim([0 2.5]);
    
%% %Calculating the angle and distance to the objects
[meh, ~] = size(clustersCentroids);    
for i = 1:meh
angle(i) = (atan2(clustersCentroids(i,2),clustersCentroids(i,1))*180)/pi
dist(i) = sqrt(clustersCentroids(i,1)^2+clustersCentroids(i,2)^2)
end
end

