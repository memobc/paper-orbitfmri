%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script runs first level background connectivity analyses in CONN
% (https://sites.google.com/view/conn/) Whitfield-Gabrieli & Nieto-Castanon (2012)

% Computes ROI-to-ROI connectivity during selected task events while
% regressing out all task- and performance-related activity.

% Analyses are used to verify the PM/AT structure and to explore changes
% in background connectivity of these networks between encoding and retrieval.

% *** 'help conn_batch' for CONN batch code options ***

% Rose Cooper - last updated Jan 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specify directories and file locations

clearvars; clc;

b.scriptdir = pwd;
addpath(b.scriptdir);

TR    = 1.5; %TR for all functional scans
scans = 466; %number of TRs per scan run

modSpace  = 'MNI-unsmoothed';
model     = 'Task-Connectivity_weighted'; %name for this connectivity project
framework = 'Features-SuccessPrecision-noEmot'; %first level univariate model to use as framework

%%%%%% define conn analysis details %%%%%
analysis_name = [model '_ROI-to-ROI'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create output directory
b.outDir = ['/gsfs0/data/ritcheym/data/fmri/orbit/analysis/orbit/connectivity/conn-models/' modSpace '/'];
if ~exist(b.outDir,'dir')
    mkdir(b.outDir);
end

% SPM/CONN toolboxes
b.spmDir = '/gsfs0/data/ritcheym/repos/toolboxes-fmri/'; %lab fmri toolboxes
addpath(genpath(b.spmDir));

% structural scans?
b.strcDir  = '/gsfs0/data/ritcheym/data/fmri/orbit/data/derivs/fmriprep/';
% functional preprocessed data? (note. smoothed folder contains copies of both smoothed and unsmoothed .niis)
b.funcDir  = '/gsfs0/data/ritcheym/data/fmri/orbit/data/derivs/smoothed/';
% ROI files?
b.ROIdir  = '/gsfs0/data/ritcheym/data/fmri/orbit/analysis/orbit/myROIs/MNI-space/';
% first level covariates (created from SPM.xX.X for conn) --> to regress out task/memory signal
b.covDir  = ['/gsfs0/data/ritcheym/data/fmri/orbit/analysis/orbit/connectivity/conn-first-level-covars/' modSpace '/' framework '/'];
% univariate first level model? (for other SPM.mat info)
b.modDir  = ['/gsfs0/data/ritcheym/data/fmri/orbit/analysis/orbit/univariate/first-level/' modSpace '/' framework '/'];

% define ROI names (PM network, Hipp, AT network):
myROIs = {'ANG','PREC','PCC','RSC','PHC','pHIPP','aHIPP','PRC','AMYG','FUS','ITC','OFC'};


clear batch
%% 1. SETUP SUBJECTS

batch.filename=fullfile(b.outDir,[analysis_name '.mat']); % New conn_*.mat project name

%get subjects who have first level covariate folders for CONN:
subjects = dir(b.covDir);
subjs = {};
for s = 1:length(subjects)
    if ~isempty(strfind(subjects(s).name,'sub')) && subjects(s).isdir
        subjs = [subjs; {subjects(s).name}];
    end
end
NSUBS = size(subjs,1);

batch.Setup.isnew     = 1; %0 if want to update existing project
batch.Setup.nsubjects = NSUBS;
batch.Setup.RT        = TR; %TR (seconds)


%% 2. GET FUNCTIONALS

fprintf('\nGetting functional scans...\n');

batch.Setup.functionals = repmat({{}},[NSUBS,1]); % initialize main functional volumes for each subject/session

% grab all niis in preprocessed functional folder, as only valid runs were copied there
for nsub=1:NSUBS
    fprintf('...%s\n',subjs{nsub});
    
    % a) main 4D functionals :
    func_regexp = ['\<' subjs{nsub} '.*task-Memory.*MNI.*_preproc.nii']; %MNI unsmoothed (no 'smooth' prefix), for ROI analyses
    clear funcFiles
    funcFiles = cellstr(spm_select('FPList', [b.funcDir subjs{nsub} '/'], func_regexp)); %find nii file for every session
    if isempty(funcFiles)
        error('no functional files found.');
    end
    nsessions = length(funcFiles);
    
    for nses=1:nsessions
        batch.Setup.functionals{nsub}{nses}{1} = funcFiles{nses,:};
    end
    
    batch.Setup.roifunctionals.roiextract = 1; % 1 indicates to use the main functional files % 4 indicates separate user defined data files (e.g. smoothed/unsmoothed)
end %note: each subject's data is defined by one single (4D) file per session


%% 3. GET STRUCTURALS

fprintf('\nGetting structural scans...\n');

batch.Setup.structurals = repmat({{}},[NSUBS,1]); %inirialize variable for structural scan names

% note. conn manual suggests can use brainmask.nii if data preprocessed using freesurfer and analyzing in native space
for nsub=1:NSUBS
    fprintf('...%s\n',subjs{nsub});
    
    % find MNI space T1 preprocessed image from fmriprep
    strc_regexp = ['\<' subjs{nsub} '.*MNI.*_preproc.nii'];
    strcFile    = cellstr(spm_select('FPList', [b.strcDir subjs{nsub} '/anat/'], strc_regexp));
    if isempty(strcFile)
        error('no structural file found.');
    end
    
    if size(strcFile,1) > 1 %already unzipped (so zipped and unzipped versions found)
        remove = cellstrfind(strcFile,'.gz');
        strcFile(remove,:) = []; %remove zipped version
        strcFile = char(strcFile);
    elseif size(strcFile,1) == 1 % need to unzip if only 1 file
        strcFile = char(strcFile);
        [~,~,ext] = fileparts(strcFile);
        % unzip if not done so already:
        if strcmp(ext,'.gz')
            gunzip(strcFile);
            strcFile = strcFile(:,1:end-length(ext)); %remove extension from target filename
        end
    end
    batch.Setup.structurals{nsub} = strcFile;
    
    %get freesurfer or fmriprep grey/WM/csf masks:
    %conn_importaseg([b.strcDir subjs{nsub} '/mri/aseg.mgz']); %use this function to get freesurfer masks --> saved as c1_aseg (grey), c2_aseg (white) and c3_aseg (csf)
    %NOTE - fmriprep output is probabilistic atlases --> CONN default is p
    %>= .5 if not binary masks, but can edit this, and can also ask CONN to mask analyses by
    %GM mask if desired.
    classFiles = {};
    for m = 1:3
        if m == 1, class = 'CSF';
        elseif m == 2, class = 'GM';
        elseif m == 3; class = 'WM';
        end
        mask_regexp = ['\<' subjs{nsub} '.*MNI.*' class '_probtissue.nii'];
        maskFile   = cellstr(spm_select('FPList', [b.strcDir subjs{nsub} '/anat/'], mask_regexp));
        if isempty(maskFile)
            error('no file found for tissue type.');
        end
        
        if size(maskFile,1) > 1 %already unzipped
            remove = cellstrfind(maskFile,'.gz');
            maskFile(remove,:) = []; %remove zipped version
            maskFile = char(maskFile);
        elseif size(maskFile,1) == 1 % need to unzip if only 1 file
            maskFile = char(maskFile);
            [~,~,ext] = fileparts(maskFile);
            % unzip if not done so already:
            if strcmp(ext,'.gz')
                gunzip(maskFile);
                maskFile = maskFile(:,1:end-length(ext)); %remove extension from target filename
            end
        end
        classFiles{m} = maskFile;
    end
    % add to batch
    batch.Setup.masks.CSF.files{nsub}   = classFiles{1};
    batch.Setup.masks.Grey.files{nsub}  = classFiles{2};
    batch.Setup.masks.White.files{nsub} = classFiles{3};
end
% specify that we only want 1 dimension
batch.Setup.masks.CSF.dimensions   = 1;
batch.Setup.masks.Grey.dimensions  = 1;
batch.Setup.masks.White.dimensions = 1;


%% 4. GET ROIs

fprintf('\nGetting ROI files...\n');

% note. subject-specific ROIs can be specified (e.g. in native space) by
% using batch.Setup.rois.files{nroi}{nsub} = [subject-specific ROI file]
for nroi = 1:length(myROIs)
    batch.Setup.rois.names{nroi} = [myROIs{nroi}];
    batch.Setup.rois.files{nroi} = [b.ROIdir myROIs{nroi} '_ROI.nii']; %bilateral file
    batch.Setup.rois.dimensions{nroi} = 1; %extract single component (mean across voxels)
    batch.Setup.rois.roiextract(nroi) = 0; %1 indicates to extract from additionally specified functional data, or same as main functionals (0).
end


%% 5. LOAD TASK REGRESSORS and HRF-weighted covariates *by run*

fprintf('\nAdding task regressors...\n');

% enters target events and all task- and memory-related covariates that
% need to be regressed out of my ROI signals
for nsub=1:NSUBS
    fprintf('...%s\n',subjs{nsub});
    
    %get number of scan runs:
    nsessions = length(batch.Setup.functionals{nsub});
    %load subject's univariate first level info
    spmfile = load([b.modDir subjs{nsub} '/SPM.mat']);
    
    % estimating background connectivity for 1) encoding and 2) retrieval
    batch.Setup.conditions.names{1} = 'Encoding';
    batch.Setup.conditions.names{2} = 'Retrieval';
    
    % add onsets and durations per session for target events, as well as
    % all covariates
    for nses = 1:nsessions
        
        % 1. add task events --> encoding and remember
        for ev = 1:2
            %get onsets for this run from concatenated data
            minTime = (scans*TR)*(nses-1); maxTime = (scans*TR)*nses;
            myOns = spmfile.SPM.Sess.U(ev).ons;
            myRows = find(myOns > minTime & myOns < maxTime); myOns = myOns(myRows);
            %adjust onset times so relative to start of run not start of entire experiment
            myOns = myOns - ((scans*TR)*(nses-1));
            %add onsets
            batch.Setup.conditions.onsets{ev}{nsub}{nses} = myOns;
            %add durations:
            batch.Setup.conditions.durations{ev}{nsub}{nses} = spmfile.SPM.Sess.U(ev).dur(myRows);
        end
        
        
        % 2. add HRF-convolved (1 x nscans) covariates for all nuisance and
        % task-related activity.
        % Note. these vectors were extracted from first level SPM.xX.X -
        % which contains HRF-convolved regressors for all events and
        % corresponding parametric modulators - and saved as individual
        % text files.
        
        % grab all covariates associated with this session for this subject:
        sessFiles = sort(cellstr(spm_select('FPList', [b.covDir subjs{nsub} '/'], ['.*_run0' num2str(nses) '.txt'])));
        
        % only need covariates for encoding and remember events, and motion
        encoding = cellstrfind(sessFiles,'/Encoding')';
        retrieval = cellstrfind(sessFiles,'/Remember')';
        motion = cellstrfind(sessFiles,'/motion')';
        sessFiles = sessFiles([encoding,retrieval,motion],:);
        
        % add names and covariate files to conn_batch
        regNames = {}; %store for denoising
        for cov = 1:length(sessFiles)
            %get regressor name
            curReg = strsplit(sessFiles{cov},'/');
            curReg = curReg{end}; %grab final part, which is the regressor name
            curReg = strsplit(curReg,'_'); %now just remove the runID from name
            curReg = curReg{1};
            regNames{cov} = curReg;
            
            %add to batch
            batch.Setup.covariates.names{cov} = curReg;
            batch.Setup.covariates.files{cov}{nsub}{nses} = sessFiles{cov};
        end
    end %end of loop through sessions
end


%% 6. DEFINE ANALYSIS steps

% CONN Setup
batch.Setup.analyses=1; %[1,2] %ROI-to-ROI(1) and seed-to-voxel(2)
batch.Setup.voxelmask=1; % 1. Explicit mask, 2. implicit subject specific
batch.Setup.voxelmaskfile=[b.ROIdir 'wb_graymatter_mask.nii']; %for any whole-brain, mask with gray matter
batch.Setup.voxelresolution=3; %same as functional volumes

% CONN Denoising --> detrending and regress out specified regressors from functional data
batch.Denoising.filter=[0.008 inf]; %matches spm high pass filtering --> recommended https://www.nitrc.org/forum/message.php?msg_id=11196   %[0.01,0.1]; %band-pass filter (Hz) - resting state default
batch.Denoising.detrending=1;       %1: linear detrending
batch.Denoising.confounds.names=regNames; %all of my task-related and nuisance covariates
for r = 1:length(regNames)
    batch.Denoising.confounds.deriv{r}=0; %do not add derviatives
end

% CONN First level
batch.Analysis.analysis_number='taskbackground'; %numerical index or custom string for analysis ID
batch.Analysis.measure=1;     % 1=bivariate correlation, 3=bivariate regression, 4=multivariate regression, see: https://www.nitrc.org/forum/message.php?msg_id=14924
batch.Analysis.weight=2;      % 1.no within-condition weighting, 2. HRF (if estimating connectivity for event-related design)
batch.Analysis.modulation=0;  % 0 = standard weighted, 1 = gPPI of condition-specific temporal modulation factor.
batch.Analysis.type=1; %1='ROI-to-ROI', 2=Seed-to-voxel, 3=all
batch.Analysis.sources=batch.Setup.rois.names; % define seeds (all source rois)


%% Run first level analyses

fprintf('\n ----- Running Setup, Denoising, and First level analyses -----\n');
batch.Setup.done          = 1; % 0 to specify files files but not run step
batch.Setup.overwrite     = 1; % 0 if updating existing project and want to keep previous files
batch.Denoising.done      = 1;
batch.Denoising.overwrite = 1;
batch.Analysis.done       = 1;
batch.Analysis.overwrite  = 1;

conn_batch(batch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%