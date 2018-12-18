function orbit_fMRI_Task(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to run Orbit fMRI Memory Task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rose Cooper - Created December 2017
%
% ENCODING:
% - participants encode a series of 'events' comprising a perspective of a
% panorama scene and an object in a specific color within the scene. 
% - also hear a negative or neutral noise with each object-scene display

% RETRIEVAL:
% - participants are show a grayscale object that they previously saw and
% are asked to remember all the features associated with the object and
% then to decide if the object is a negative (bomb) or neutral (safe) object.
% - participants are then asked to re-create the object color and
% background scene location.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs = filesep;

%folders for memory task stimuli:
scenePath = [S.stimPath 'scene_views' fs];
objPath   = [S.stimPath 'objects' fs];
soundPath = [S.stimPath 'sounds' fs];

%create data file name for this subject:
S.FileName = sprintf('sub%s_%s_startRun%d_Orbit.mat',S.subNum,S.Date,S.startRun);


%% Define stimuli numbers

%task structure details:
stimIterations   = 120; %360 degree circle for responses but actually jumps 3 degrees when key pressed (for speed), so 120 unique angles.
stimChange       = 360/stimIterations; % change around circular spaces per key press, in degrees
TotalNumTrials   = 144; %number of 'events' (unique scene-object pairings)
numBlocks        = 6;   %number of blocks for study-test.
TotalTrialsBlock = TotalNumTrials/numBlocks;
numPanos         = 6;   %number of unique scene panoramas
numSounds        = 12;  %number of unique sounds
% get start and end trial per block:
blkFirst         = (1:numBlocks)*TotalTrialsBlock - (TotalTrialsBlock - 1);
blkEnd           = (1:numBlocks)*TotalTrialsBlock;


try % -------------------------------------------------------------------
%% Load in stimuli and specify timings:

textColor  = 255; %white

% load in pre-generated stimuli presentation orders
load([S.taskPath 'Orders_fMRI.mat'],'StudyOrd','TestOrd');
% get the orders specific to this subject (as indexed by S.PresOrd), and
% remove headers
S.studyOrder = StudyOrd{S.PresOrd}(2:end,:);
S.testOrder  = TestOrd{S.PresOrd}(2:end,:);

% read in list of panorama and sound names
[~,sceneNames,~] = xlsread([S.stimPath 'Scene_List.xlsx']);
[~,soundNames,~] = xlsread([S.stimPath 'Sound_List.xlsx']);

%display times:
if S.debug ==0
debug_factor = 1;
else %if debugging
debug_factor = 20; %task speeds X times faster if in debugging/testing mode 
end

getReadyTime  = 12/debug_factor;  %between study and test phase and for feedback
fixTime       = 1/debug_factor;   %fixation ITI between trials
encTime       = 6/debug_factor;   %encoding time for item-context display
rememberTime  = 4/debug_factor;   %time to remember everything during the test phase
emotionTime   = 2/debug_factor;   %time to respond to emotion question
respTime      = 6/debug_factor;   %time to provide memory judgements for color and scene


%load in image for inter-trial fixation:
fixFile  = [S.taskPath 'fix.jpg'];
fix      = imread(fixFile);
fix      = imresize(fix,[S.screenH*0.1 S.screenH*0.1]);   
trialfix = Screen('MakeTexture',S.window,fix);
fixSize  = size(fix);
fixRect  = [S.xCenter-fixSize(2)/2 S.yCenter-fixSize(1)/2 S.xCenter+fixSize(2)/2 S.yCenter+fixSize(1)/2];


% load objects in study trial order (all blocks)-----------------:
for ev = 1:TotalNumTrials 
    message = strcat('Loading objects. \n\nPlease wait...',sprintf('working on %d/%d',ev,TotalNumTrials));
    DrawFormattedText(S.window,message,'center','center',textColor);
    Screen('Flip', S.window);  %'flips' message to the screen
    
    im_file = strcat(objPath,S.studyOrder{ev,4}); %grab object file name
  
    %save LAB version for later color manipulation
    [the_image,~,alp]  = imread(im_file); 
    the_image    = imresize(the_image,[240 240]);
    alpha{ev}    = imresize(alp,[240 240]);
    savedLab{ev}(:,:,:) = colorspace('rgb->lab', the_image); % convert to lab

    %save grayscale version for remember question - note texture saved in
    %study order *not* test order
    gray_image = rgb2gray(the_image);
    gray_image = imresize(gray_image,[300 300]);
    alphaGray  = imresize(alp,[300 300]);
    gray_image(:,:,2) = alphaGray;
    grayObjects{ev} = gray_image;  
   
end %end of loop through objects
% define rectangles for object display
objSize = size(the_image);
objRect = [S.xCenter-objSize(2)/2 S.yCenter-objSize(1)/2 S.xCenter+objSize(2)/2 S.yCenter+objSize(1)/2];
graySize = size(gray_image);
grayRect = [S.xCenter-graySize(2)/2 S.yCenter-graySize(1)/2 S.xCenter+graySize(2)/2 S.yCenter+graySize(1)/2];


% load scenes - in scene list order (*not* trial unique) ------------:
myScenes = {};
for ev = 1:numPanos
  message = strcat('Loading scenes. \n\nPlease wait...',sprintf('working on %d/%d',ev,numPanos));
  DrawFormattedText(S.window,message,'center','center',textColor);
  Screen('Flip', S.window);  %'flips' message to the screen

  im_dir = strcat(scenePath,sprintf('%s.jpg',sceneNames{ev,1}),fs); %get name of scene panorama
  
 % read in all unique 120 images (perspectives) for scene, but don't make individual textures:
 for i = 1:stimIterations
  im_file = strcat(im_dir,sprintf('%d.jpeg',i));  
  the_image = imread(im_file);
  the_image = imresize(the_image,[600 800]);  
  the_image = addborder(the_image, 3, [255 255 255], 'inner'); 
  myScenes(i,ev) = {the_image};
 end %end of loop through perspectives
    
end %end of loop through scene contexts
imSize = size(the_image); %all images are same size so just uses last one
imRect = [S.xCenter-imSize(2)/2 S.yCenter-imSize(1)/2 S.xCenter+imSize(2)/2 S.yCenter+imSize(1)/2];


% load sounds - in sound list order (*not* trial unique) -------------------:
freq       = 44100;
nrchannels = 1;
SHandle    = PsychPortAudio('Open',[],[],1,freq,nrchannels);   % opens audio player
buffer = [];

for ev = 1:numSounds
  message = strcat('Loading sounds. \n\nPlease wait...',sprintf('working on %d/%d',ev,numSounds));
  DrawFormattedText(S.window,message,'center','center',textColor);
  Screen('Flip', S.window);  %'flips' message to the screen
  
  soundFile = strcat(soundPath,soundNames{ev,1},'.wav');  %get name of sound file
  [Swavedata, infreq]  = audioread(soundFile);            %load .wav file
  
  %resample freq if different from target
  if infreq ~= freq
  Swavedata = resample(Swavedata, freq, infreq);       
  end
  %create buffer for each sound:
  buffer(ev) = PsychPortAudio('CreateBuffer', [], Swavedata'); 
end


% define location for object within scene:
X = S.xCenter;
Y = S.yCenter + (imSize(1)/2) - (objSize(1)/2) - (imSize(1)/10);
RectLocation = CenterRectOnPointd(objRect, X, Y);

%location of object for remember test question
RememberObjLocation = CenterRectOnPointd(grayRect, X, (S.yCenter-200));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through study-test blocks - START FMRI RUN

for blk = S.startRun:S.endRun %just in case need to restart or run subset of blocks
%% ----------------------- RUN ENCODING --------------------------------%

message = sprintf('STUDY PHASE - BLOCK %d \n\nPlease press ''s'' to start.',blk);
DrawFormattedText(S.window,message,'center','center',textColor);
Screen('Flip',S.window);

% wait for key press to start experiment - press BEFORE start scan sequence
while 1
   [keyDown, ~, keyCodes] = KbCheck(-1);
   if keyDown
       if keyCodes(S.startKey)  %wait for experimenter to trigger
           break;
       elseif keyCodes(S.quitKey)
           sca; return;
       end
   end
end 
Priority(MaxPriority(S.window));


%now set KbQueueCheck ready to receive scanner pulse trigger
triglist = zeros(1,256);     %%create a list
triglist(S.triggerKey) = 1;  %%set target keys to 1
KbQueueCreate(S.trigDevice,triglist); %%make cue


 % -------- WAIT FOR SCANNER PULSE ---------- %
message = 'Get Ready...';
DrawFormattedText(S.window,message,'center','center',textColor);
Screen('Flip',S.window);
S.startTimeWait(blk) = GetSecs; % time when started waiting for scanner trigger

if S.scanner == 0 %if behavioral-only just wait a few seconds before start

WaitSecs(getReadyTime/2);

elseif S.scanner == 1 %if running in scanner...
    
KbQueueStart(S.trigDevice);   %start listening for trigger key!    
KbQueueWait(S.trigDevice);    %wait for trigger key (will only respond to trigger due to KbQueueCreate set up)
KbQueueRelease(S.trigDevice); %when trigger found, clear queue object

end
S.startTimeEnc(blk) = GetSecs; % < ----- ! start of experiment timing - triggered by scanner
goTime = S.startTimeEnc(blk); %counter to make sure all stim displays stay to schedule


if S.scanner == 1
%%% Restrict key input to subject's button box keys (0-4):
RestrictKeysForKbCheck([S.leftKey, S.rightKey, S.threeKey, S.fourKey, S.responseKey])
end


for ev = blkFirst(blk):blkEnd(blk)  % < -------- start of loop through encoding events

    goTime = goTime + fixTime;  
    %%% Fixation before encoding trial
    Screen('DrawTexture', S.window, trialfix,[],fixRect);
    onset_fixEnc = Screen('Flip',S.window);
    %store relative to block start time:
    S.onset_fixEnc(ev) = onset_fixEnc - S.startTimeEnc(blk);
    
% ---------------------------------------------------------------------- %
%build stim loading into fix time so computation doesn't add time:
    
    %define scene and sound for this trial based on study order
    curScene  = S.studyOrder{ev,5}; %1-X identity of scene panorama for this trial
    curSound  = S.studyOrder{ev,7}; %1-X identity of sound for this trial
    PsychPortAudio('FillBuffer', SHandle, buffer(curSound)); %fill sound buffer for this trial

    %create color object:
    colAng   = S.studyOrder{ev,13}; %color angle to display
    newImage = RotateImage(savedLab{ev}, (colAng*stimChange));% set object to this color
    newImage(:,:,4) = alpha{ev};
    curobjTexture   = Screen('MakeTexture', S.window, newImage);
    
    %create scene perspective:
    sceAng  = S.studyOrder{ev,14}; %scene angle
    cursceTexture = Screen('MakeTexture', S.window, myScenes{sceAng,curScene});  
    
    %.....SCENE ..... +
    Screen('DrawTexture', S.window, cursceTexture, [], imRect);      
    %.... COLOR OBJECT +
    Screen('DrawTexture', S.window, curobjTexture, [], RectLocation);   

% ---------------------------------------------------------------------- %
while GetSecs < goTime  % ------ > fixed time
WaitSecs(0.001);
end

    % duration of encoding event:
    goTime = goTime + encTime;  
       
    %.....SOUND....... play 1 repetition
    PsychPortAudio('Start', SHandle, 1);
    
    %show object + scene
    onset_objEnc = Screen('Flip',S.window);
    %store onset relative to block start time:
    S.onset_objEnc(ev) = onset_objEnc - S.startTimeEnc(blk);
while GetSecs < goTime  % ------ > fixed time
WaitSecs(0.001);
end

    %close textures:
    Screen('Close', curobjTexture);
    Screen('Close', cursceTexture);
    

 %----------------------------------------------------------------%  
end %end of loop through events

    goTime = goTime + fixTime; 
    %%% Fixation after all encoding trials
    Screen('DrawTexture', S.window, trialfix,[],fixRect);
    onset_fixEndEnc = Screen('Flip',S.window);
    S.onset_fixEndEnc(blk) = onset_fixEndEnc - S.startTimeEnc(blk);
    
    %save/update subject's data file:
    save([S.sDir S.FileName],'S');
    
while GetSecs < goTime  % ------ > fixed time
WaitSecs(0.001);
end

  
    
%% ----------------------- RUN RETRIEVAL --------------------------------%

goTime = goTime + getReadyTime; 
message = sprintf('TEST PHASE - BLOCK %d \n\nGet ready to start...',blk);
DrawFormattedText(S.window,message,'center','center',textColor);
startTimeRet = Screen('Flip',S.window);
S.startTimeRet(blk) = startTimeRet - S.startTimeEnc(blk); % < ----- start of retrieval timing relative to run start
while GetSecs < goTime  % ------ > fixed time
WaitSecs(0.001);
end


for ev = blkFirst(blk):blkEnd(blk)  % -------- > start of loop through test events

goTime = goTime + fixTime; 
%%% fixation before remember question:
Screen('DrawTexture', S.window, trialfix,[],fixRect);
onset_fixRemember = Screen('Flip',S.window);
%store onset relative to start time:
S.onset_fixRemember(ev) = onset_fixRemember - S.startTimeEnc(blk);

%define matching study trial number (as objects pre-stored in this order)
curTrial  = S.testOrder{ev,1};

%make texture for grayscale object:
grayTexture = Screen('MakeTexture', S.window, grayObjects{curTrial}); %gray version of object

while GetSecs < goTime  % ------ > fixed time
WaitSecs(0.001);
end


% 1a. REMEMBER ALL ------------------------------------------------------

goTime = goTime + rememberTime; 

%display cue to remember as much as possible about the object
Screen('DrawTexture', S.window, grayTexture, [], RememberObjLocation);
message = 'Remember...';
DrawFormattedText(S.window,message,'center',S.yCenter + 100,textColor);

onset_remember = Screen('Flip',S.window);
%store onset relative to start time:
S.onset_remember(ev) = onset_remember - S.startTimeEnc(blk);

% ---------------------------------------------------------------------- %
%%% Make object and scene images during this remember time to speed up
%%% subsequent rotation

    %define features for this test trial based on randomized test order
    qOrd    = cell2mat(S.testOrder(ev,11:12));%counterbalanced Q order (color then scene; scene then color)
    colAng  = S.testOrder{ev,15};   %color start test angle
    sceAng  = S.testOrder{ev,16};   %scene start test angle
    
    %define scene for this test trial based on randomized test order
    curScene  = S.testOrder{ev,5}; %1-X identity of scene panorama
    %make image texture for first test perspective:
    cursceTexture = Screen('MakeTexture', S.window, myScenes{sceAng,curScene});  
    
    %load in all color versions of this object:
    myObjects = {};
    savedLab1   = savedLab{curTrial}(:,:,:);
    alpha1      = alpha{curTrial};    
for c = 1:stimIterations
    newImage        = RotateImage(savedLab1, (c*stimChange));
    newImage(:,:,4) = alpha1(:,:);    
    myObjects(c)    = {newImage}; 
end
   
    %make image texture for start test color:
    curobjTexture  = Screen('MakeTexture', S.window, myObjects{colAng}); 

%%% ---------------------------------------------------------------------
while GetSecs < goTime  % ------ > fixed time
WaitSecs(0.001);
end


% 1b. EMOTION: ---------------------------------------------------------

goTime = goTime + emotionTime; 

%display emotion screen:
Screen('DrawTexture', S.window, grayTexture, [], RememberObjLocation);
message = 'Bomb or Safe?';
DrawFormattedText(S.window,message,'center',S.yCenter + 100,textColor);
message = '1 = Sure BOMB     2 = Maybe BOMB     3 = Maybe SAFE     4 = Sure SAFE';
DrawFormattedText(S.window,message,'center',S.yCenter + 250,textColor);
    
onset_emotion = Screen('Flip',S.window);
%store onset relative to start time:
S.onset_emotion(ev) = onset_emotion - S.startTimeEnc(blk);

RT = -1; key = -1;
while GetSecs < goTime
    [keyDown, pressTime, keyCode] = KbCheck(-1); %waits for a key press
    if keyDown && RT == -1 %first response counts
        if keyCode(S.leftKey), key = 1; RT = pressTime - onset_emotion;
        elseif keyCode(S.rightKey), key = 2; RT = pressTime - onset_emotion;
        elseif keyCode(S.threeKey), key = 3; RT = pressTime - onset_emotion;
        elseif keyCode(S.fourKey), key = 4; RT = pressTime - onset_emotion;
        end
    end
    WaitSecs(0.001); %prevent overload
end

S.emotionRT(ev)  = RT;
S.emotionKey(ev) = key;

Screen('Close', grayTexture);


%. 2 PRECISION QUESTIONS: ----------------------------------------------

for q = 1:2 %loop through both questions
 
 if qOrd(q)==1 %if 1 first, ask color, if not, ask scene first (counterbalanced)
     
 %%% 2A. COLOR RETRIEVAL %%%%%%:
 
 goTime = goTime + fixTime; 
 %fixation
 Screen('DrawTexture', S.window, trialfix,[],fixRect);
 onset_fixCol = Screen('Flip',S.window);
 S.onset_fixCol(ev) = onset_fixCol - S.startTimeEnc(blk);
while GetSecs < goTime  % ------ > fixed time
WaitSecs(0.001);
end


goTime = goTime + respTime; 

DrawFormattedText(S.window,'Color','center',(S.yCenter + (imSize(1)/2) + 50),textColor);
%draw images:
Screen('DrawTexture', S.window, cursceTexture, [], imRect);
Screen('DrawTexture', S.window, curobjTexture, [], RectLocation);
onset_col  = Screen('Flip',S.window);
%convert to relative to start time to store:
S.onset_col(ev) = onset_col - S.startTimeEnc(blk);

curAngle = colAng;
RT = -1;
while GetSecs < goTime  % ------ > fixed time
    
    oldAng = curAngle; %previous angle before this iteration:
    [keyIsDown,~,keyCodes] = KbCheck(-1);  %checks if key pressed and which one
    if keyIsDown
        if keyCodes(S.rightKey) && RT == -1 % only move if response not confirmed yet
            if curAngle < stimIterations %increase in angle if moving right
                curAngle = curAngle + 1;
            else, curAngle = 1;  %reset to 1 if moved beyond 360
            end
        elseif keyCodes(S.leftKey) && RT == -1 % only move if response not confirmed yet
            if curAngle > 1 %decrease in angle if moving left
                curAngle = curAngle - 1;
            else, curAngle = stimIterations;  %reset to 360 if moved below 1
            end
        elseif keyCodes(S.responseKey) && RT == -1 % confirmed response -- first RT counts
            RT = GetSecs - onset_col;
        end
    end
             
      %change color if moved:
      if curAngle ~= oldAng
      DrawFormattedText(S.window,'Color','center',(S.yCenter + (imSize(1)/2) + 50),textColor);
     
      %draw scene image according to start test or chosen angle (if already done question):
      Screen('DrawTexture', S.window, cursceTexture, [], imRect);
      
      %change object color based on current angle:
      Screen('Close', curobjTexture);
      curobjTexture = Screen('MakeTexture', S.window, myObjects{curAngle}); 
      Screen('DrawTexture', S.window, curobjTexture, [], RectLocation);
      Screen('Flip',S.window);
      end
      
      WaitSecs(0.001); %prevent overload
        
end % End of while loop---------------

   S.ColResp(ev) = curAngle;
   S.ColRT(ev)   = RT;


elseif qOrd(q)==2

  %%% 2B. SCENE RETRIEVAL: %%%%%%%%
  
  goTime = goTime + fixTime; 
  %fixation
  Screen('DrawTexture', S.window, trialfix,[],fixRect);
  onset_fixSce = Screen('Flip',S.window);
  S.onset_fixSce(ev) = onset_fixSce - S.startTimeEnc(blk);
while GetSecs < goTime  % ------ > fixed time
WaitSecs(0.001);
end


goTime = goTime + respTime; 

DrawFormattedText(S.window,'Scene','center',(S.yCenter + (imSize(1)/2) + 50),textColor);
%draw images:
Screen('DrawTexture', S.window, cursceTexture, [], imRect);
Screen('DrawTexture', S.window, curobjTexture, [], RectLocation);
onset_sce = Screen('Flip',S.window);  
%store onset relative to start time:
S.onset_sce(ev) = onset_sce - S.startTimeEnc(blk);

curAngle = sceAng;
RT = -1;
while GetSecs < goTime  % ------ > fixed time
    
    oldAng = curAngle; %previous angle before this iteration:
    [keyIsDown,~,keyCodes] = KbCheck(-1);  %checks if key pressed and which one
    if keyIsDown
        if keyCodes(S.rightKey) && RT == -1 % only move if response not confirmed yet
            if curAngle < stimIterations  %increase in angle if moving right
                curAngle = curAngle + 1;
            else, curAngle = 1; %reset to 1 if moved beyond 360
            end
        elseif keyCodes(S.leftKey) && RT == -1 % only move if response not confirmed yet
            if curAngle > 1 %decrease in angle if moving left
                curAngle = curAngle - 1; 
            else, curAngle = stimIterations; %reset to 360 if moved below 1
            end
        elseif keyCodes(S.responseKey) && RT == -1 % confirmed response -- first RT counts
            RT = GetSecs - onset_sce;
        end
    end
            
      %change scene perspective if moved:
      if curAngle ~= oldAng
      DrawFormattedText(S.window,'Scene','center',(S.yCenter + (imSize(1)/2) + 50),textColor);
      
      %change scene perspective:
      Screen('Close', cursceTexture);
      cursceTexture = Screen('MakeTexture', S.window, myScenes{curAngle,curScene}); 
      Screen('DrawTexture', S.window, cursceTexture, [], imRect);

      %draw object in chosen/start test color:
      Screen('DrawTexture', S.window, curobjTexture, [], RectLocation);
      Screen('Flip',S.window);
      end
      
      WaitSecs(0.001); %prevent overload
      
end % End of while loop---------------

   S.SceResp(ev) = curAngle;
   S.SceRT(ev)   = RT;


 end % end of 'if' statement for counterbalance order

end % end of loop through precision questions
 
    Screen('Close', curobjTexture);
    Screen('Close', cursceTexture);
 
    
end %end of loop through events -----------------------------------------

    goTime = goTime + fixTime; 
    %%% Fixation before feedback
    Screen('DrawTexture', S.window, trialfix,[],fixRect);
    onset_fixFeedback = Screen('Flip',S.window);
    S.onset_fixFeedback(blk) = onset_fixFeedback - S.startTimeEnc(blk); 
while GetSecs < goTime  % ------ > fixed time
WaitSecs(0.001);
end


%% present feedback on block just completed
    
    goTime = goTime + getReadyTime; 

    % 1. Emotion memory for this block - percentage correct
    for t = blkFirst(blk):blkEnd(blk)
        if S.testOrder{t,7} <=6 %negative sound
        if S.emotionKey(t) ==1 || S.emotionKey(t) ==2 %if said a bomb
            emotCorr(t) = 1;
        else, emotCorr(t) = 0;
        end
        elseif S.testOrder{t,7} >=7 %if a neutral sound
        if S.emotionKey(t) ==3 || S.emotionKey(t) ==4 %if said safe
            emotCorr(t) = 1;
        else, emotCorr(t) = 0;
        end
        end
    end
    emoPercent = round(mean(emotCorr(blkFirst(blk):blkEnd(blk)))*100);
    emoMessage = sprintf('You correctly identified bombs and safe objects on %d%% of trials.',emoPercent);
    DrawFormattedText(S.window,emoMessage,'center',S.yCenter - 150,[200 0 50]);
    
    
    % 2. Color memory (good if <=45 degrees - just to give idea of 'accuracy'):
    studyAngles = cell2mat(S.testOrder(blkFirst(blk):blkEnd(blk),13)); %encoded angles
    respAngles  = S.ColResp(blkFirst(blk):blkEnd(blk))'; %responses
    
    %convert to 360 and then to radians to calculate error:
    studyAngles = wrap((studyAngles*stimChange)/180*pi);
    respAngles  = wrap((respAngles*stimChange)/180*pi);  
    resultsCol  = abs(wrap(respAngles - studyAngles)); %absolute difference between response and study
    resultsCol  = round(rad2deg(resultsCol));
    
    colPercent = round((length(find(resultsCol <= 45))/length(resultsCol))*100);
    colMessage = sprintf('You were close to the object''s original color on %d%% of trials.',colPercent);
    DrawFormattedText(S.window,colMessage,'center',S.yCenter - 50,[0 200 50]);
        

    % 3. Scene memory (good if <=45 degrees - just to give idea of 'accuracy'):
    studyAngles = cell2mat(S.testOrder(blkFirst(blk):blkEnd(blk),14)); %encoded angles
    respAngles  = S.SceResp(blkFirst(blk):blkEnd(blk))'; %responses
    
    %convert to 360 and then to radians to calculate error:
    studyAngles = wrap((studyAngles*stimChange)/180*pi);
    respAngles  = wrap((respAngles*stimChange)/180*pi);  
    resultsSce  = abs(wrap(respAngles - studyAngles)); %absolute difference between response and study
    resultsSce  = round(rad2deg(resultsSce));
    
    scePercent = round((length(find(resultsSce <= 45))/length(resultsSce))*100);
    sceMessage = sprintf('You were close to the object''s original scene location on %d%% of trials.',scePercent);
    DrawFormattedText(S.window,sceMessage,'center',S.yCenter + 50,[0 50 200]);
    
    if blk < numBlocks %don't need to give additional feedback on last block
    %custom phrase based on performance:
      if emoPercent >= 80 && colPercent >= 70 && scePercent >= 70 
          message = 'Great! Keep up the good work.';
      elseif emoPercent >= 70 && colPercent >= 50 && scePercent >= 50 
          message = 'Good job! Please try to remember the object details as accurately as possible.';
      else
          message = 'Please try to remember the object details as accurately as possible.';
      end
    DrawFormattedText(S.window,message,'center',S.yCenter + 200,textColor);
    end
    

    %flip feedback to screen
    onset_Feedback = Screen('Flip',S.window);
    S.onset_Feedback(blk) = onset_Feedback - S.startTimeEnc(blk);
    
    %save data at end of each retrieval block:
    save([S.sDir S.FileName],'S');

while GetSecs < goTime  % ------ > fixed time
WaitSecs(0.001);
end


% --------------------------------------------------------------- %
% ----------------- SCAN SEQUENCE STOPS HERE -------------------- %
% --------------------------------------------------------------- %

    
%%% Restore all key input:
RestrictKeysForKbCheck([])
%%%

    
end %end of loop through blocks

catch ME  %if an error, close psychtoolbox
Screen('Close',S.window);
fprintf('Error: %s\n',ME.message)
end %end of try loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%% ADDITIONAL FUNCTIONS
% ----------------------------------------------------------------------%
% by Tim Brady - Rotates color of objects in LAB space
% ----------------------------------------------------------------------%
function newRgb = RotateImage(lab, r)
  x = lab(:,:,2);
  y = lab(:,:,3);
  v = [x(:)'; y(:)'];
  vo = [cosd(r) -sind(r); sind(r) cosd(r)] * v;
  lab(:,:,2) = reshape(vo(1,:), size(lab,1), size(lab,2));
  lab(:,:,3) = reshape(vo(2,:), size(lab,1), size(lab,2));
  newRgb = colorspace('lab->rgb', lab) .* 255;
end

% ----------------------------------------------------------------------%
% by Paul Bays - wraps angles in circular space
% ----------------------------------------------------------------------%
function X = wrap(Y, bound)
% X = WRAP(Y)
%   Maps Y values onto circular space (-PI <= X < PI). 
%   X = WRAP(Y, BOUND) specifies alternative bounds (default is PI).
%   --> www.paulbays.com
if nargin<2, bound = pi; end
X = mod(Y + bound, bound*2) - bound;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%