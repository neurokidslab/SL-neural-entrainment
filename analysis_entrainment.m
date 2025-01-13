clear, close all
restoredefaultpath
clc

%% ------------------------------------------------------------------------
%% ADD TOOLBOXES

% EEGLAB
pp = pwd;
cd(fullfile('path_to_eeglab','eeglab2022.1'))
eeglab
close all
cd(pp);

% APICE
addpath(genpath(fullfile('path_to_apice','eeg_preprocessing')))

% functions for the entrainmet
addpath(genpath(fullfile('path_to_entrainment_codes', 'codes')))

%% ------------------------------------------------------------------------
%% PARAMETERS

% Filter
% -------------------------------------------------------------------------
filt_highpass = 0.2;

% Define BT before segmenting
% -------------------------------------------------------------------------
BCall.nbt           = [0.70 0.50 0.30];  % define bat channels
BTall.nbc           = [0.70 0.50 0.30];  % define bad times
BTall.maskTime      = 0.250;        % also mark as bad surronding samples it was 0.5
BTall.minGoodTime   = 2.000;        % short intervals between bad segments will be marked as bad
BTall.minBadTime    = 0.100;        % too short periods will not be considered as bad

% Short epoch
% -------------------------------------------------------------------------
Epoch.event = {'eRNl' 'eSTl' 'eSTs'};  % conditions to epoch
Epoch.nwords = 1;  % short epoch numper of words
Epoch.tw = [0, 0.50];  % short epoch time windows
epochlenght = 0.50;

% Short epochs: Define Bad Segments and Electrodes
BCall_s.nbt           = [0.30 0.30];    % define bad channels during all recording
BCep_s.nbt            = [0.30 0.10/epochlenght];   % define bad channels per epoch
BTep_s.nbc            = [0.30 0.30];    % define bat times
BTep_s.minGoodTime    = 1.00;           % shorter intervals between bad segments will be marked as bad
BTep_s.minBadTime     = 0.10;           % too short periods will not be considered as bad
BTep_s.maskTime       = 0.00;           % also mark as bad surronding samples

% Short epochs: Define bad epoch
% based on the amount of bad data
DefBEa_s.maxloops     = 1;
DefBEa_s.limBCTa      = 1.00;              % maximun bad data per epoch
DefBEa_s.limBTa       = 0.10/epochlenght;  % maximun bad times per epoch
DefBEa_s.limBCa       = 0.30;   % maximun channles with some bad data 
DefBEa_s.limCCTa      = 1.00;   % maximun interpolated data
% based on distance to the average
DefBEdT_s.limDist     = 2.0;
DefBEdT_s.limBadDist  = 0.10/epochlenght; 
DefBEdT_s.maxloops    = 1;
DefBEdT_s.rmvmean     = 1;

% Long epochs
% -------------------------------------------------------------------------
nWords = 15;
epochlenght = nWords*0.50;

% Long epochs: Define Bad Segments and Electrodes
BCall_l.nbt           = [0.30 0.30];             % define bad channels during all recording
BCep_l.nbt            = [0.30 0.10/epochlenght]; % define bad channels per epoch
BTep_l.nbc            = [0.30 0.30];             % define bad times
BTep_l.minGoodTime    = 1.00;                    % shorter intervals between bad segments will be marked as bad
BTep_l.minBadTime     = 0.10;                    % too short periods will not be considered as bad
BTep_l.maskTime       = 0.00;                    % also mark as bad surronding samples

% Long epochs: Define bad epochs
% based on the amount of bad data
DefBEa_l.maxloops     = 1;
DefBEa_l.limBCTa      = 1.00;  % maximun bad data per epoch
DefBEa_l.limBTa       = 0.10/epochlenght;   % maximun bad times per epoch 
DefBEa_l.limBCa       = 0.30;   % maximun channles with some bad data 
DefBEa_l.limCCTa      = 1.00;   % maximun interpolated data
% based on distance to the average
DefBEdT_l.limDist     = 2.0;
DefBEdT_l.limBadDist  = 0.10/epochlenght; 
DefBEdT_l.maxloops    = 1;
DefBEdT_l.rmvmean     = 1;

% To reshape the short epoch into log epochs after rejection
reshape.f 		= {'eventtype' 'eventRep'};
reshape.fkeep   = {'eventtype','eventBlock','eventBlock60','eventCnd','eventRep'};
reshape.nWords  = nWords;
reshape.tlims   = [0 500];

% Artifacts correction
% -------------------------------------------------------------------------
Int.Spl.p = 0.3;
Int.Spl.pneigh  = 1;

% Artifatcs detection 
% -------------------------------------------------------------------------
ArtMot_s = [];
ArtMot_s(1).algorithm        = 'eega_tRejMaxChangeEpoch';
ArtMot_s(1).loops            = [1];
ArtMot_s(1).P.refdata        = 1;
ArtMot_s(1).P.refbaddata     = 'replacebynan';
ArtMot_s(1).P.dozscore       = 0;
ArtMot_s(1).P.thresh         = 2.0;
ArtMot_s(1).P.relative       = 1;
ArtMot_s(1).P.xelectrode     = 0;

ArtMot_l = [];
ArtMot_l(1).algorithm        = 'eega_tRejMaxChangeEpoch';
ArtMot_l(1).loops            = [1];
ArtMot_l(1).P.refdata        = 1;
ArtMot_l(1).P.refbaddata     = 'replacebynan';
ArtMot_l(1).P.dozscore       = 0;
ArtMot_l(1).P.thresh         = 2.0;
ArtMot_l(1).P.relative       = 1;
ArtMot_l(1).P.xelectrode     = 0;

ArtPwr = [];
ArtPwr(1).algorithm        = 'eega_tRejPwr';
ArtPwr(1).loops            = [1];
ArtPwr(1).P.refdata        = 0;
ArtPwr(1).P.dozscore       = 1;
ArtPwr(1).P.twdur          = nWords*0.50;
ArtPwr(1).P.twstep         = 1;
ArtPwr(1).P.frqband        = [1 10];    %(Hz) frequency bands
ArtPwr(1).P.thresh         = [-Inf 2]; % threshold for the power 
ArtPwr(1).P.relative       = [1]; 

% Entrainment
% -------------------------------------------------------------------------

% Define what to remove (epochs artifacts)
exclude_bc = 1;  % exclude channels marked as bad

% conditions to put toghether
factorscnd      = {'eventCnd'};     
factorstyp      = {'eventtype'};   
ordcnd          = [1 2];          % [2 1 3] | [4 3 1 2]
ordtype         = [1 2 3];          % [2 1 3] | [4 3 1 2]

% to estimate the power
P.foutput         = [1/(0.50*nWords) 8]; % output frequnecies

% SNR computation
ftarget = [1.33 2 4];
freqSNR = 0.8;
noutSNR = 2.0;
typeSNR = 'zscore'; 
fbandoutput = [];

% time windows analysis
nepxtimew = 16;  % epochs per time window
nepstep = 2;     % time step

%% ------------------------------------------------------------------------
%% RUN THE ANALYSIS


% Load the data
% -------------------------------------------------------------------------
folder_data_in = 'data_example';
file_data_in = 'a_i_is_a_fh1e-01_fl40_Session 20180926 S01.set';
EEG = pop_loadset(file_data_in,folder_data_in);

% Add events to segment the data
% -------------------------------------------------------------------------
% rename the first DIN1 trigger after the first RNl as eRNl
EEG = eega_changeevents(EEG, 'DIN1',  'type', 'eRNl', 'RNl', 1, 'first', [0  Inf]);
% rename the first DIN1 trigger after the second RNl as eRNl
EEG = eega_changeevents(EEG, 'DIN1',  'type', 'eRNl', 'RNl', 2, 'first', [0  Inf]);
% rename the first DIN1 trigger after STl as eSTl
EEG = eega_changeevents(EEG, 'DIN1',  'type', 'eSTl', 'STl', 1, 'first', [0  Inf]);
% rename the first DIN1 trigger after each STs as eSTs
EEG = eega_changeevents(EEG, 'DIN1',  'type', 'eSTs', 'STs', 1, 'first', [0  Inf]);
EEG = eega_changeevents(EEG, 'DIN1',  'type', 'eSTs', 'STs', 2, 'first', [0  Inf]);
EEG = eega_changeevents(EEG, 'DIN1',  'type', 'eSTs', 'STs', 3, 'first', [0  Inf]);
EEG = eega_changeevents(EEG, 'DIN1',  'type', 'eSTs', 'STs', 4, 'first', [0  Inf]);
EEG = eega_changeevents(EEG, 'DIN1',  'type', 'eSTs', 'STs', 5, 'first', [0  Inf]);
EEG = eega_changeevents(EEG, 'DIN1',  'type', 'eSTs', 'STs', 6, 'first', [0  Inf]);
EEG = eega_changeevents(EEG, 'DIN1',  'type', 'eSTs', 'STs', 7, 'first', [0  Inf]);
EEG = eega_changeevents(EEG, 'DIN1',  'type', 'eSTs', 'STs', 8, 'first', [0  Inf]);
EEG = eega_changeevents(EEG, 'DIN1',  'type', 'eSTs', 'STs', 9, 'first', [0  Inf]);
EEG = eega_changeevents(EEG, 'DIN1',  'type', 'eSTs', 'STs', 10, 'first', [0  Inf]);
% remove all DIN1 events
EEG = eega_removeevent(EEG, 'DIN1', [], 'all');
% add events for each word (every 0.5 s)
EEG = eega_addmultipleevents(EEG, 'eRNl', 'eRNl', (Epoch.nwords*0.5:Epoch.nwords*0.5:119.5));
EEG = eega_addmultipleevents(EEG, 'eSTl', 'eSTl', (Epoch.nwords*0.5:Epoch.nwords*0.5:119.5));
EEG = eega_addmultipleevents(EEG, 'eSTs', 'eSTs', (Epoch.nwords*0.5:Epoch.nwords*0.5:29.5));
% add description for different events
EEG = eega_changeevents(EEG, 'eRNl' , 'Cnd', 'RND', 'RNl', 1, (1:240/Epoch.nwords),     [0  Inf] ); 
EEG = eega_changeevents(EEG, 'eSTl' , 'Cnd', 'STR', 'STl', 1, (1:240/Epoch.nwords),     [0  Inf] );
EEG = eega_changeevents(EEG, 'eSTs' , 'Cnd', 'STR', 'STs', 1, (1:600/Epoch.nwords),     [0  Inf] );
EEG = eega_changeevents(EEG, 'eRNl', 'Rep', 'rep1', 'RNl', 1, (1:240/Epoch.nwords),     [0 121] );
EEG = eega_changeevents(EEG, 'eSTl', 'Rep', 'rep1', 'STl', 1, (1:240/Epoch.nwords),     [0 121] );
EEG = eega_changeevents(EEG, 'eSTs', 'Rep', 'rep2', 'STs', 1, (1:(60*10)/Epoch.nwords),  [0  Inf] );
EEG = eega_changeevents(EEG, 'eRNl', 'Block', 'b1', 'RNl', 1, (1:240/Epoch.nwords),     [0 121] );
EEG = eega_changeevents(EEG, 'eSTl', 'Block', 'b1', 'STl', 1, (1:240/Epoch.nwords),     [0 121] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block', 'b1', 'STs', 1, (1:60/Epoch.nwords),      [0  31] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block', 'b2', 'STs', 2, (1:60/Epoch.nwords),      [0  31] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block', 'b3', 'STs', 3, (1:60/Epoch.nwords),      [0  31] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block', 'b4', 'STs', 4, (1:60/Epoch.nwords),      [0  31] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block', 'b5', 'STs', 5, (1:60/Epoch.nwords),      [0  31] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block', 'b6', 'STs', 6, (1:60/Epoch.nwords),      [0  31] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block', 'b7', 'STs', 7, (1:60/Epoch.nwords),      [0  31] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block', 'b8', 'STs', 8, (1:60/Epoch.nwords),      [0  31] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block', 'b9', 'STs', 9, (1:60/Epoch.nwords),      [0  31] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block', 'b10', 'STs', 10, (1:60/Epoch.nwords),     [0  31] );
EEG = eega_changeevents(EEG, 'eRNl', 'Block60', 'b1', 'RNl', 1, 'all',  [-0.1 59.9] );
EEG = eega_changeevents(EEG, 'eRNl', 'Block60', 'b2', 'RNl', 1, 'all',  [59.9 119.9] );
EEG = eega_changeevents(EEG, 'eSTl', 'Block60', 'b1', 'STl', 1, 'all',  [-0.1 59.9] );
EEG = eega_changeevents(EEG, 'eSTl', 'Block60', 'b2', 'STl', 1, 'all',  [59.9 119.9] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block60', 'b1', 'STs', 1, 'all',  [-0.1  29.9] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block60', 'b1', 'STs', 2, 'all',  [-0.1  29.9] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block60', 'b2', 'STs', 3, 'all',  [-0.1  29.9] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block60', 'b2', 'STs', 4, 'all',  [-0.1  29.9] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block60', 'b3', 'STs', 5, 'all',  [-0.1  29.9] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block60', 'b3', 'STs', 6, 'all',  [-0.1  29.9] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block60', 'b4', 'STs', 7, 'all',  [-0.1  29.9] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block60', 'b4', 'STs', 8, 'all',  [-0.1  29.9] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block60', 'b5', 'STs', 9, 'all',  [-0.1  29.9] );
EEG = eega_changeevents(EEG, 'eSTs', 'Block60', 'b5', 'STs', 10, 'all', [-0.1  29.9] );
% remove unecessary events
EEG = eega_removeevent(EEG, 'RNl',  [], 'all');
EEG = eega_removeevent(EEG, 'STl', [], 'all');
EEG = eega_removeevent(EEG, 'STs', [], 'all');

% Short epochs
% -------------------------------------------------------------------------
% high pass filter
EEG = pop_eegfiltnew(EEG, filt_highpass, [], [], 0, [], [], 0);
% define bad channels and times on continuos data
EEG = eega_tDefBTBC(EEG, BTall.nbc,BCall.nbt,BCall.nbt,'keeppre',0,'minBadTime',BTall.minBadTime,'minGoodTime',BTall.minGoodTime,'maskTime',BTall.maskTime);
% epoch
EEG = eega_epoch(EEG, Epoch.event, Epoch.tw);
EEG = eega_definefactors(EEG, {'eventtype' 'eventBlock' 'eventBlock60' 'eventCnd' 'eventRep'});


% Deal with artefacts
% -------------------------------------------------------------------------
% further detect big signal chnages in short epochs
EEG = eega_tArtifacts(EEG, ArtMot_s, 'KeepRejPre', 1,'FilterLp',[],'FilterHp',[]);
% define bad times and channnles per short epoch
EEG = eega_tDefBTBC(EEG, BTep_s.nbc,BCep_s.nbt,BCall_s.nbt,'keeppre','BT','minBadTime',BTep_s.minBadTime,'minGoodTime',BTep_s.minGoodTime,'maskTime',BTep_s.maskTime);
% define bad short epochs
EEG = eega_tDefBEbaddata(EEG, DefBEa_s,'keeppre',0,'plot', 0);
EEG = eega_tDefBEdist(EEG, DefBEdT_s.limDist, DefBEdT_s.limBadDist, 'maxloops', DefBEdT_s.maxloops,'keeppre',1,'rmvmean',DefBEdT_s.rmvmean, 'plot', 0);


% Entrainment per contidion
% -------------------------------------------------------------------------

% remove artefacted hort epochs
eeg = eega_rmvbadepochs(EEG);
% interpolate bad channels per epoch
eeg = eega_tInterpSpatialEEG(eeg, Int.Spl.p,'pneigh', Int.Spl.pneigh);
% re-detect remaining big changes in the signal
eeg = eega_tArtifacts(eeg, ArtMot_s, 'KeepRejPre', 1,'FilterLp',[],'FilterHp',[]);
% reshape into log epochs
eeg = dft_reshape(eeg, reshape.nWords, reshape.f, reshape.fkeep, reshape.tlims);
% detect big signal chnages in long epochs
eeg = eega_tArtifacts(eeg, ArtPwr, 'KeepRejPre', 1,'FilterLp',[],'FilterHp',[]);
% define bad times and channnles per long epoch
eeg = eega_tDefBTBC(eeg, BTep_l.nbc,BCep_l.nbt,BCall_l.nbt,'keeppre','BT','minBadTime',BTep_l.minBadTime,'minGoodTime',BTep_l.minGoodTime,'maskTime',BTep_l.maskTime);
% define bad long epochs
eeg = eega_tDefBEbaddata(eeg, DefBEa_l,'keeppre',0,'plot', 0);
eeg = eega_tDefBEdist(eeg, DefBEdT_l.limDist, DefBEdT_l.limBadDist,'maxloops', DefBEdT_l.maxloops,'keeppre',1,'rmvmean',DefBEdT_l.rmvmean, 'plot', 0);
% remove bad long epoch
eeg = eega_rmvbadepochs(eeg);
% everage reference
eeg = eega_refavg(eeg);
% data normalization
eeg = eega_normalization(eeg, 'epochs','single','electrodes','all','mean',0,'sd',[],'BadData','replacebynan');

% compute neural entrainment per condition
ent = dft_applyent(eeg, factorscnd, P, exclude_bc);

% SNR computation    
freqs = ent.freqs;
nbinsnorm_1 = 0;
nbinsnorm_2 = round( freqSNR / (freqs(2)-freqs(1)) );
binsSNR = [(-nbinsnorm_2:-(1+nbinsnorm_1)) ((1+nbinsnorm_1):nbinsnorm_2)];
ent.plvsnr = dft_norm( ent.plv, ent.freqs,...
    'binsnorm', binsSNR,'fband',fbandoutput,'ftarget',ftarget,'typeSNR',typeSNR,'noutliers',noutSNR);


% Entrainment time course
% -------------------------------------------------------------------------

% interpolate bad channels per epoch
eeg = eega_tInterpSpatialEEG(EEG, Int.Spl.p,'pneigh', Int.Spl.pneigh);
% re-detect remaining big changes in the signal
eeg = eega_tArtifacts(eeg, ArtMot_s, 'KeepRejPre', 1,'FilterLp',[],'FilterHp',[]);

% compute neural entrainment along time
ent_tc = dft_applyenttimecourse(eeg, nepxtimew, nepstep,...
                                'reshape',  reshape,....
                                'exclude_bc', exclude_bc,...
                                'BCall',    BCall_l,...
                                'BCep',     BCep_l,...
                                'BTep',     BTep_l,...
                                'DefBEa',   DefBEa_l,...
                                'DefBEdT',  DefBEdT_l,...
                                'ent',      P);

% SNR computation
freqs = ent.freqs;
nbinsnorm_1 = 0;
nbinsnorm_2 = round( freqSNR / (freqs(2)-freqs(1)) );
binsSNR = [(-nbinsnorm_2:-(1+nbinsnorm_1)) ((1+nbinsnorm_1):nbinsnorm_2)];
ent_tc.plvsnr = dft_norm( ent_tc.plv, ent_tc.freqs,...
    'binsnorm', binsSNR,'fband',fbandoutput,'ftarget',ftarget,'typeSNR',typeSNR,'noutliers',noutSNR);

