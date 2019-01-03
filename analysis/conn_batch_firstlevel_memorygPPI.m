%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Batch script for running gPPI analyses using the CONN toolbox
% (https://sites.google.com/view/conn/) Whitfield-Gabrieli & Nieto-Castanon (2012)

% Requires that each subject already has a specified univariate first level
% model (SPM.mat).
% From this, task/memory-related covariates were saved as text files
% for CONN gPPI: the toolbox cannot extract parametric modulators from an SPM.mat,
% only event regressors, so parametric modulators have to be entered separately as
% 1xnscan covariates using the HRF-convolved vectors in SPM.xX.X.

% PIPELINE:
% 1. Specify subjects
% 2. Load functional (unsmoothed and/or smoothed) .niis --> smoothed for
% whole brain (if applicable) and unsmoothed for ROI analyses per subject
% and session.
% 3. Load structural scan and GM/WM/CSF masks per subject.
% 4. Add ROI nii files.
% 5. Define task events if applicable.
% 6. Define first level covariates - load parametric modulators and other nuisance
% regressors as text files from SPM.xX.X.
% 7. Specify denoising and analysis steps
% 8. Run first level analysis

% *** 'help conn_batch' for CONN batch code options ***

% Rose Cooper - last updated Jan 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specify directories and file locations

clearvars; clc;

b.scriptdir = pwd;
addpath(b.scriptdir);

TR    = 1.5; %TR for all functional scans
scans = 466; %number of TRs per scan run


% define the first level GLM to extract pmod information from for current CONN model
modSpace = 'MNI-unsmoothed';
model    = 'Features-SuccessPrecision-noEmot';

%%%%%% define conn analysis details %%%%%
analysis_name = [model '_ROI-to-ROI'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create output directory
b.outDir = ['/gsfs0/data/ritcheym/data/fmri/orbit/analysis/orbit/connectivity/conn-models/' modSpace '/'];
if ~exist(b.outDir,'dir')
    mkdir(b.outDir);
end

% SPM/CONN toolboxes
% **NOTE - the lab SPM 12 repo has implicit threshold changed from standard 0.8 default to -inf default (so analyses all voxels at first level -- include ROI or whole-brain masks)
b.spmDir = '/gsfs0/data/ritcheym/repos/toolboxes-fmri/'; %lab fmri toolboxes
addpath(genpath(b.spmDir));

% structural scans?
b.strcDir  = '/gsfs0/data/ritcheym/data/fmri/orbit/data/derivs/fmriprep/';
% functional preprocessed data? (note. smoothed folder contains both smoothed and unsmoothed copies of niis)
b.funcDir  = '/gsfs0/data/ritcheym/data/fmri/orbit/data/derivs/smoothed/';
% ROI files?
b.ROIdir  = '/gsfs0/data/ritcheym/data/fmri/orbit/analysis/orbit/myROIs/MNI-space/';
% univariate first level models? (for SPM.mat info)
b.modDir  = ['/gsfs0/data/ritcheym/data/fmri/orbit/analysis/orbit/univariate/first-level/' modSpace '/' model '/'];
% first level covariates (created for conn)?
b.covDir  = ['/gsfs0/data/ritcheym/data/fmri/orbit/analysis/orbit/connectivity/conn-first-level-covars/' modSpace '/' model '/'];

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

batch.Setup.isnew     = 1; %1;  %0 if want to update existing analysis
batch.Setup.nsubjects = NSUBS;
batch.Setup.RT        = TR; %TR (seconds)


%% 2. GET FUNCTIONALS

fprintf('\nGetting functional scans...\n');

batch.Setup.functionals = repmat({{}},[NSUBS,1]); %main (smoothed) functional volumes for each subject/session
batch.Setup.roifunctionals.roiextract_functionals = repmat({{}},[NSUBS,1]); %additional (unsmoothed) functional volumes for ROI analyses

% grab all niis in preprocessed functional folder, as only valid runs were copied there
for nsub=1:NSUBS
    fprintf('...%s\n',subjs{nsub});
    
    % a) main functionals (for whole-brain if running):
    func_regexp = ['\<smooth.*task-Memory.*MNI.*_preproc.nii']; %smoothed
    clear funcFiles
    funcFiles = cellstr(spm_select('FPList', [b.funcDir subjs{nsub} '/'], func_regexp)); %find subjects's nii file for every session
    if isempty(funcFiles)
        error('no smoothed functional files found.');
    end
    nsessionsSm = length(funcFiles);
    
    % add to batch
    for nses=1:nsessionsSm
        batch.Setup.functionals{nsub}{nses}{1} = funcFiles{nses,:};
    end
    
    % b) additional for calculating ROI means:
    func_regexp = ['\<' subjs{nsub} '.*task-Memory.*MNI.*_preproc.nii']; %MNI unsmoothed (no 'smooth' prefix), for ROI analyses
    clear funcFiles
    funcFiles = cellstr(spm_select('FPList', [b.funcDir subjs{nsub} '/'], func_regexp)); %find nii file for every session
    if isempty(funcFiles)
        error('no unsmoothed functional files found.');
    end
    nsessionsUn = length(funcFiles);
    
    % add to batch
    for nses=1:nsessionsUn
        batch.Setup.roifunctionals.roiextract_functionals{nsub}{nses}{1} = funcFiles{nses,:};
    end
    batch.Setup.roifunctionals.roiextract = 1; % 1 indicates same as main functional files % 4 indicate user defined data files
    
    if nsessionsSm ~= nsessionsUn
        error('Different number of smoothed and unsmoothed data files found.');
    end
end %note: each subject's data is defined by one single (4d) file per session


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


%% 5. LOAD TASK REGRESSORS *by run*

fprintf('\nAdding task regressors...\n');

% get all of this information from each subject's univariate first level
% SPM.mat file:
% Note. if just analyzing variables that are pmods, need to add 'dummy'
% events (covering entire scan) to modulate with memory covariates for gPPI.
% see: https://www.nitrc.org/forum/message.php?msg_id=14924
for nsub=1:NSUBS
    fprintf('...%s\n',subjs{nsub});
    
    %load subject's univariate first level info
    spmfile = load([b.modDir subjs{nsub} '/SPM.mat']);
    
    % Note. only looking at 'remember' events.
    regNames = {};
    nEvents = length(spmfile.SPM.Sess.U);
    for ev = 2
        regNames = [regNames, spmfile.SPM.Sess.U(ev).name];
    end
    
    %get number of scan runs:
    nsessions = length(batch.Setup.functionals{nsub});
    
    nameCount = 0; % count through univariate regressor names
    regCount  = 0; % count number of gPPI regressors
    %add onsets and durations for events
    for ev = 2 %only 'Remember'
        nameCount = nameCount + 1;
        % get number of parametric modulators associated with this event
        numP = length(spmfile.SPM.Sess.U(ev).name)-1;
        
        if numP > 0 %only interested in modeling modulator memory effects if pmods present
            for p = 1:numP
                nameCount = nameCount + 1;
                regCount = regCount + 1;
                batch.Setup.conditions.names{regCount} = regNames{nameCount};
                batch.Setup.conditions.param(regCount) = regCount;  % indicate the covariate index to modulate this 'event'
                
                % add event onsets and durations per session
                for nses = 1:nsessions
                    %event onsets - start of session
                    batch.Setup.conditions.onsets{regCount}{nsub}{nses} = 0;
                    %event durations - whole session as using session length
                    %covariate modulator
                    batch.Setup.conditions.durations{regCount}{nsub}{nses} = inf;
                end
            end
        end
    end
end


%% 6. ADD HRF-weighted COVARIATES/PMODS to modulate above event blocks

% NOTE - the condition 'weights' used by conn are the original parametric
% modulators values - CONN does not change or mean-center, so make sure the
% raw values are mean-centered appropriately (here across all runs) already.
fprintf('\nAdding first level covariates...\n');

%already generated from create_txtcovars_conn
for nsub=1:NSUBS
    fprintf('...%s\n',subjs{nsub});
    
    spmfile = load([b.modDir subjs{nsub} '/SPM.mat']);
    nEvents = length(spmfile.SPM.Sess.U); %number of event regressors
    
    %get number of scan runs:
    nsessions = length(batch.Setup.functionals{nsub});
    for nses = 1:nsessions
        regCount = 0;
        %add onsets and durations for events
        for ev = 2
            if length(spmfile.SPM.Sess.U(ev).name) > 1 %if this event has a pmod, covariate is parametric modulator
                for p = 2:length(spmfile.SPM.Sess.U(ev).name) %start from 2 to avoid first main event regressor
                    regCount = regCount + 1;
                    pname = spmfile.SPM.Sess.U(ev).name{p};
                    batch.Setup.covariates.names{regCount} = pname;
                    % get text file for this covariate, for this subject and run
                    batch.Setup.covariates.files{regCount}{nsub}{nses} = [b.covDir subjs{nsub} '/' pname '_run0' num2str(nses) '.txt'];
                end
            end
        end
        %add covariates for denoising (FD, aCompCorr, and 6 movement parameters):
        regCount = regCount + 1;
        batch.Setup.covariates.names{regCount} = 'denoising';
        batch.Setup.covariates.files{regCount}{nsub}{nses} = [b.covDir subjs{nsub} '/motion_run0' num2str(nses) '.txt'];
    end %end of loop through sessions
end %end of loop through subjects


%% 7. DEFINE ANALYSIS steps
% CONN Setup
batch.Setup.analyses=[1,2]; %ROI-to-ROI(1) and seed-to-voxel(2)
batch.Setup.voxelmask=1; % 1. Explicit mask, 2. implicit subject specific
batch.Setup.voxelmaskfile=[b.ROIdir 'wb_graymatter_mask.nii'];
batch.Setup.voxelresolution=3; %same as functional volumes
batch.Setup.outputfiles=[0,0,1,0,0,0]; %saves out seed to voxel r maps (see conn_batch for further details) for whole-brain analyses

% CONN Denoising --> detrending and regress out nuisance regressors 'motion' from preprocessed functional data
batch.Denoising.filter=[0.008 inf]; %matches spm high pass filtering --> recommended https://www.nitrc.org/forum/message.php?msg_id=11196 %[0.01,0.1]; % band-pass filter (in Hz) - resting state default
batch.Denoising.detrending=1;      % 1: linear detrending
batch.Denoising.confounds.names={'denoising'};  %no need to remove task effects here as gPPI controls for these in regression model
batch.Denoising.confounds.deriv{1}=0;  %do not add derviatives
% ^^ Note. further adding WM and CSF initiates aCompCorr method (recommended),
% whereas adding GM performs global signal regression (not recommended). I
% have already included aCompCorr as a regressor in 'denoising' generated
% by fMRIprep, so don't need to further run via CONN.
% CONN also automatically de-means the functional data of each run, to
% remove differences in average signal between runs (as analyses are always
% done on conctenated data) --> see 'results/preprocessing/'

% CONN First level
batch.Analysis.analysis_number='gppi'; %numerical index or custom string for analysis ID
batch.Analysis.measure=3;   %1=bivariate correlation, 3=bivariate regression, 4=multivariate regression, see: https://www.nitrc.org/forum/message.php?msg_id=14924
batch.Analysis.weight=1;    %1=none, 2=HRF % don't need HRF weighting as not modeling events -- pmods already HRF-convolved.
batch.Analysis.modulation=1; % 0 = standard weighted, 1 = gPPI of condition-specific temporal modulation factor.
batch.Analysis.type=3; %1='ROI-to-ROI', 2=Seed-to-voxel, 3=all
batch.Analysis.sources=batch.Setup.rois.names; % define seeds (all source rois)
% CONN defaults to entering all conditions into the same analysis (full gPPI model)


%% 8. Run first level analyses

fprintf('\n ----- Running Setup, Denoising, and First level analyses -----\n');
batch.Setup.done          = 1; % 0 to specify files files but not run step
batch.Setup.overwrite     = 1; % 0 if updating existing project and want to keep previous files
batch.Denoising.done      = 1;
batch.Denoising.overwrite = 1;
batch.Analysis.done       = 1;
batch.Analysis.overwrite  = 1;

conn_batch(batch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%