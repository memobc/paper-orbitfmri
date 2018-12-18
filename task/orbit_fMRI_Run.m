function orbit_fMRI_Run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rose Cooper - Created December 2017

%%% Master script to run the Orbit fMRI experiment %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc;

fs = filesep;
myPath = [pwd fs];

S.taskPath  = [myPath 'task' fs 'fMRI' fs]; %task scripts
S.stimPath  = [myPath 'stimuli' fs 'fMRI' fs]; %stimuli
S.dataPath  = [myPath 'data' fs 'fMRI' fs]; %to save data

rng('shuffle');


%% Enter Participant/Session details:

prompt = {'Date (MMDDYY):', ...
    'Subject number:',...
    'Start run:',...
    'End run:',...
    'Scanner? (1 = fMRI, 0 = behav-only)',...
    'Debugging? (default=0):'};
defaults = {'XXXX18','0XX','X','1','6','1','0'};
answer = inputdlg(prompt,'Experimental setup',1,defaults);

[S.Date, S.subNum, S.startRun, S.endRun, S.scanner, S.debug] = deal(answer{:});

S.startRun = str2num(S.startRun);
S.endRun   = str2num(S.endRun);
S.scanner  = str2num(S.scanner);
S.debug    = str2num(S.debug);

%define unique stimuli combinations/trial list based on subject number
S.PresOrd  = str2num(S.subNum);

%create directory for subject data:
S.sDir	= [S.dataPath S.subNum fs];
mkdir(S.sDir);



%% Specify Psychtoolbox set up:

%basic set up and configures the keyboard using KbName('UnifyKeyNames'):
PsychDefaultSetup(1);
Screen('Preference', 'SkipSyncTests', 1);

%initialization of the sound driver:
PsychPortAudio('Close');
InitializePsychSound(1);

%screen defaults:
screenColor   = 0;   %black background
textColor     = 255; %white text
defaultText   = 38;  %font size

HideCursor; %hides cursor from screen

% Get screenNumber to display stimuli (should be 0 if behavioral-only and 1 if scanner and detects a second screen).
screens   = Screen('Screens');
screenNum = max(screens);

% open full-screen psychtoolbox window:
[S.window, winRect] = PsychImaging('OpenWindow', screenNum, screenColor, []);

%set text properties:
Screen('TextFont', S.window, 'Arial');
Screen('TextSize', S.window, defaultText);

%get center coordinates and dimensions of screen (used for displaying stimuli):
[S.xCenter, S.yCenter] = RectCenter(winRect);
S.screenW = winRect(3); S.screenH = winRect(4);

Screen('BlendFunction', S.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%%% Make sure all keys initialized at start:
RestrictKeysForKbCheck([])
%%%



%% Specify response devices/keys:

%keys for experimenter to start and quit:
S.startKey    = KbName('s');    % 's' to start each block
S.quitKey     = KbName('q');    % 'q' to quit at at block

%define response keys for participant:
if S.scanner == 0 %if behavioral
    %behavioral keys
    S.leftKey     = KbName('LeftArrow');
    S.rightKey    = KbName('RightArrow');
    S.responseKey = KbName('space');
    
elseif S.scanner == 1 %if fMRI -- use button box
    S.triggerKey  = KbName('=+'); % scanner pulse is '=' --> used to trigger experiment
    %button box keys are 0-4
    S.leftKey     = KbName('1!');
    S.rightKey    = KbName('2@');
    S.threeKey    = KbName('3#');
    S.fourKey     = KbName('4$');
    S.responseKey = KbName('0)');
    
    %get scanner trigger device ID:
    deviceString='Celeritas Dev'; %fibre optic scanner response box
    [id,name] = GetKeyboardIndices;
    
    S.trigDevice=0;
    for i=1:length(name) %for each possible device
        if strcmp(name{i},deviceString) %compare the name to the name you want
            S.trigDevice=id(i); %grab the correct id, and exit loop
            break;
        end
    end
    if S.trigDevice==0 %%error if can't find button box/trigger device
        error('No device by that name was detected');
    end
end %end of loop to define response keys



%% RUN EXPERIMENT

% 1. Run Orbit Memory task: ---------------------------------------------
orbit_fMRI_Task(S)

message = 'End of Experiment!';
DrawFormattedText(S.window,message,'center','center',textColor);
Screen('Flip',S.window);
WaitSecs(1);
KbStrokeWait(-1);
Screen('Close',S.window);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%