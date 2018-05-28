%This sript is made by ROB18_Gr465 AAU 

%This scipt is a modified version of the loop_back sript provided by
%play_rec. The Original loop_back sript recorded the data on the microphone
%and play the recorded sample on the speaker.
%Inorder to run this sript it is required that playrec is install.

%This script is made in such a way that the speaker emits a puls of sound and 
%the microphones records the pulse and the echo coming back 

function [recBuffer] = Play_rec_new3( playDeviceID, recDeviceID, playChanList, recChan, Fs)
%LOOP_BACK Loops back an input channel onto output channel(s)

%% setup
pageSize = 1024;    %size of each page processed
time = 0.1;         %record time
v_sound = 343.2;    %speed of sound

% Setup microphone positions (x-axis displacement)

mic_1x = -0.15;     %Position of mic 1 on x-axis
mic_2x = 0;         %Position of mic 2 on x-axis
mic_3x = 0.15;      %Position of mic 3 on x-axis

runMaxSpeed = true;  %When true, the processor is used much more heavily 
                     %(ie always at maximum), but the chance of skipping is 
                     %reduced

%% 
% Check if the number of play Channels is real row vector  
if ~isreal(playChanList) || length(playChanList) < 1 ...
    || ndims(playChanList)~=2 || size(playChanList, 1)~=1

    error ('playChanList must be a real row vector with at least 1 element');
end

%Test if current initialisation is ok
if(playrec('isInitialised'))
    if(playrec('getSampleRate')~=Fs)
        fprintf('Changing playrec sample rate from %d to %d\n', playrec('getSampleRate'), Fs);
        playrec('reset');
    elseif(playrec('getPlayDevice')~=playDeviceID)
        fprintf('Changing playrec play device from %d to %d\n', playrec('getPlayDevice'), playDeviceID);
        playrec('reset');
    elseif(playrec('getRecDevice')~=recDeviceID)
        fprintf('Changing playrec record device from %d to %d\n', playrec('getRecDevice'), recDeviceID);
        playrec('reset');       
    elseif(playrec('getPlayMaxChannel')<max(playChanList))
        fprintf('Resetting playrec to configure device to use more output channels\n');
        playrec('reset');
    elseif(playrec('getRecMaxChannel')<max(recChan))
        fprintf('Resetting playrec to configure device to use more input channels\n');
        playrec('reset');
    end
end

%Initialize if not initialized
if(~playrec('isInitialised'))
    fprintf('Initialising playrec to use sample rate: %d, playDeviceID: %d and recDeviceID: %d\n', Fs, playDeviceID, recDeviceID);
    playrec('init', Fs, playDeviceID, recDeviceID, max(playChanList), max(recChan))
    
    pause(0.1);
end
    
if(~playrec('isInitialised'))
    error ('Unable to initialise playrec correctly');
end

if(playrec('pause'))
    fprintf('Playrec was paused - clearing all previous pages and unpausing.\n');
    playrec('delPage');
    playrec('pause', 0);
end

%% setup varibles 
pageNumList = [];
Flipdata = []; 
nextPageSamples = zeros(pageSize, 1);   
recBuffer = [];                         %Initialize recording array (storing of the recording)
pulse = zeros(pageSize,1);              %Initialize pulse array

%% Setting up the pulse

i = 1;
pulse(i) = 1;
pulse(i+1) = 0;
pulse(i+2) = -1;
pulse(i+3) = 0;


%% Recording

for repeatCount = 1:(time*Fs/pageSize)
    recBuffer = [recBuffer,Flipdata];
    pageNumList = [pageNumList playrec('playrec', pulse, playChanList, -1, recChan)];  
    
    if(repeatCount==1)
        % This is the first time through so reset the skipped sample count
        playrec('resetSkippedSampleCount');
        % Resetting the pulse, to only emit once in the loop
        pulse = zeros(pageSize,1);
    end
    
    if(runMaxSpeed)
        
        while(playrec('isFinished', pageNumList(1)) == 0)
        end
    else
        %Check if the page is done
        playrec('block', pageNumList(1));
    end
    %Get recordings
    nextPageSamples = playrec('getRec', pageNumList(1));
    playrec('delPage', pageNumList(1));
    Flipdata = nextPageSamples';
    pageNumList = pageNumList(2:end);
end
playrec('delPage');
    