%This Matlab script is made by ROB18 Gr465 AAU.

%Room impluse tracking simulation.

%In order to run this script it is required that RIR_generator made by
%Emanuël Habets is installed.
%https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator


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

%This Matlab script calculates the room impluse respones of a given room
%and Uses multilateration to locate the sound source inside the room.
%This is done 100 times in a for loop to simulate a tracking of the object

clear
close all
%% Obj calculator

angles = 60 * pi/180;
dists = 0.7;

%Move mic increment
inc = 0.000;

%% Setup of Variables
%-----------Reciever 1 Location-----------%
rx1 = 1.85; ry1 = 3 ; rz1 = 1 ;
%-----------Reciever 2 Location-----------%
rx2 = 2; ry2 = 3 ; rz2 = 1 ;
%-----------Reciever 3 Location-----------%
rx3 = 2.15; ry3 = 3 ; rz3 = 1 ;
%-----------Source Location-----------%
sx = rx2+(dists*cos(angles)); sy = ry2+(dists*sin(angles)); sz = 1;
%-----------Room size-----------%
rsx = 6.5 ; rsy = 5.25 ; rsz = 4 ;

%Position of mics
mic_1x = (rx2-rx1)/2;
mic_2x = 0;
mic_3x = (rx2-rx3)/2;

%% source position and angle
sourcep = [sx-rx2, sy-ry2]

%% %-----------Mandatory Inputs-----------%
c = 343.2 ; %velocity (m/s)
fs =  96000; %frequency (Hz)
s = [sx sy sz]; %Source Vector (m)
L = [rsx rsy rsz]; %Room size Vector (m)
beta = 0.3; % Reveberation time
t = 1/fs;
%% %-----------Tracking-----------%
potatodegree = [];
potatodistance = [];
%% %-----------Signal-----------%
x = zeros(1, 8000);
x(3:10)=1;
SNR=2;
%% %-----------Detail Inputs-----------%
nsample = 8000 ; %Number of samples
%mtype = 'omnidirectional'; %Microphone type
%order = -1; %Maximum reflection order (0 no rev, 1 only 1 bounce)
%dim = 3; %Room dimension (set to 2 or 3)
%orientation = [0 pi/2] ; %Direction of mic, use angle in rads
%hp_filter = ;  %Highpass filter true or false

%% Calculation and plotting

cord = [];
angle = [];
dist = [];    
counter = 0;

for q=1:100
    rx1 = rx1 + inc;  rx2 = rx2 + inc;  rx3 = rx3 + inc;
    r1 = [rx1 ry1 rz1]; %Reciever 1 location matrix (m)
    r2 = [rx2 ry2 rz2]; %Reciever 2 location matrix (m)
    r3 = [rx3 ry3 rz3]; %Reciever 3 location matrix (m)
    r = [r1;r2;r3]; %Reciever matrix (m)

%% %-----------RIR Function-----------%
h1 = rir_generator(c, fs, r1, s, L, beta, nsample);
h2 = rir_generator(c, fs, r2, s, L, beta, nsample);
h3 = rir_generator(c, fs, r3, s, L, beta, nsample);
%% %-----------Noise Generator-----------%
y1 = filter(h1,1,x);
y2 = filter(h2,1,x);
y3 = filter(h3,1,x);
Y1 = awgn(y1,SNR,'measured');
Y2 = awgn(y2,SNR,'measured');
Y3 = awgn(y3,SNR,'measured');
%% %-----------Time variable-----------%
t1 = (0:length(Y1)-1)/fs;
t2 = (0:length(Y2)-1)/fs;
t3 = (0:length(Y3)-1)/fs;

%% %Offset mics 
x12 = -3+mic_1x:0.001:3+mic_1x;
x23 = -3+mic_3x:0.001:3+mic_3x;
x13 = -3+mic_2x:0.001:3+mic_2x;

mic_1x = mic_1x-inc;
mic_2x = mic_2x-inc;
mic_3x = mic_3x-inc;

%% %Removing every thing after the intial pulse
[~,I1] = max(Y1);
[~,I2] = max(Y2);
[~,I3] = max(Y3);

for i = I1+5:length(Y1)
    Y1(i) = 0;
end
for i = I2+5:length(Y2)
    Y2(i) = 0;
end
for i = I3+5:length(Y3)
    Y3(i) = 0;
end


%% %-----------Correlation-----------%


[acor1,lag1]= xcorr(Y1,Y2);
[acor2,lag2]= xcorr(Y1,Y3);
[acor3,lag3]= xcorr(Y2,Y3);


%% %-----------Calculation of TDOAs-----------%
[~,I1] = max(abs(acor1));
lagDiff1 = lag1(I1);
timeDiff1 = lagDiff1*t;

[~,I2] = max(abs(acor2));
lagDiff2 = lag2(I2);
timeDiff2 = lagDiff2*t;

[~,I3] = max(abs(acor3));
lagDiff3 = lag3(I3);
timeDiff3 = lagDiff3*t;

%% %-----------Distances-----------%
d1=[r1;r2];
d2=[r1;r3];
d3=[r2;r3];
d1=pdist(d1,'euclidean');
d2=pdist(d2,'euclidean');
d3=pdist(d3,'euclidean');
lengtha1 = timeDiff1 * c;
lengtha2 = timeDiff2 * c;
lengtha3 = timeDiff3 * c;
%% %-----------Hyperbolas-----------%
x = -3:0.001:3;

D_12 = sqrt((rx1-rx2)^2);
D_23 = sqrt((rx2-rx3)^2);
D_13 = sqrt((rx1-rx3)^2);

hyper12 = [];
hyper23 = [];
hyper13 = [];

hyper12 = [hyper12; 1/2*sqrt(((4*x12.^2-lengtha1^2)*(D_12^2-lengtha1^2))/lengtha1^2)];
hyper23 = [hyper23; 1/2*sqrt(((4*x23.^2-lengtha3^2)*(D_23^2-lengtha3^2))/lengtha3^2)];
hyper13 = [hyper13; 1/2*sqrt(((4*x13.^2-lengtha2^2)*(D_13^2-lengtha2^2))/lengtha2^2)];

%% %intersections
%Find all intersections between all hyperbolas
[XOUT12a, YOUT12a] = intersections(x,hyper12,x,hyper23,true);
[XOUT23a, YOUT23a] = intersections(x,hyper23,x,hyper13,true);
[XOUT13a, YOUT13a] = intersections(x,hyper12,x,hyper13,true);

XOUT = [real(XOUT12a);real(XOUT23a);real(XOUT13a)];
YOUT = [real(YOUT12a);real(YOUT23a);real(YOUT13a)];

Coord = [XOUT YOUT];

%Removing all intersections with less length than 0.2m and lager than 5m
[A,~] = size(Coord);
for i = 1:(A)
   if sqrt(Coord(i,1)^2+Coord(i,2)^2) > 0.2 && sqrt(Coord(i,1)^2+Coord(i,2)^2) < 5
      newa_Coord(i,:) = Coord(i,:);
   end
end

%Removing all zeros
[row,~] = find(newa_Coord(:,1) ~= 0)
new_coord = newa_Coord(row,:)

%Quick fix for clustering
new_coord = [new_coord ; 0 , -10];
new_coord

%% %Clustering
%Find all intersection that have a length less than 0.05m to each other.
%and make a cluster that indicates a object.
[clustersCentroids, clustersGeoMedians, clustersXY] = clusterXYpoints(new_coord, 0.05, 3, 'point', 'merge');

%Storing of clusters, calculating the angle and calculating the distance.
if ~isempty(clustersCentroids)
cord = [cord; clustersCentroids]
angle = [angle; (atan2(clustersCentroids(1,2),clustersCentroids(1,1))*180)/pi]
dist = [dist; sqrt((sourcep(1)-clustersCentroids(1,1))^2+(sourcep(2)-clustersCentroids(1,2))^2)]

%Counter
counter = counter+1

%% %Plot Hyperbolas and object points
end
figure(1)
clf
plot(x,hyper12);
hold on
plot(x,hyper23);
hold on
plot(x,hyper13);
hold on
scatter(cord(:,1),cord(:,2),'filled', 'k')
hold on
scatter(sourcep(1),sourcep(2),'filled', 'g')
xlim([-1 1]);ylim([0 2]);

end