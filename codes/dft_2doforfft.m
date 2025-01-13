function [EEG] = dft_2doforfft(EEG,varargin)

%% ------------------------------------------------------------------------
%% Parameters

% Reshape
P.reshape.f         = [];
P.reshape.fkeep     = [];
P.reshape.nWords    = [];
P.reshape.tlims     = [];

% Define Bad Segments and Electrodes
P.BCBTkeeppre         = 0;
P.BCall.nbt           = [0.70 0.50 0.30 0.30];
P.BCep.nbt            = [0.70 0.50 0.30 0.10];   %0.10/diff(Epoch.tw)
P.BTep.nbc            = [0.70 0.50 0.30 0.30];
P.BTep.minGoodTime    = 1.00;             % shorter intervals between bad segments will be marked as bad
P.BTep.minBadTime     = 0.10;             % too short periods will not be considered as bad
P.BTep.maskTime       = 0.00;                % also mark as bad surronding samples

% Define bad epoch based on the amount of bad data absolute threshold
P.DefBEa.maxloops     = 1;
P.DefBEa.limBCTa      = 1.00; %0.05*(1/diff(OPT.epoch.tw));   % maximun bbad data per epoch
P.DefBEa.limBTa       = 0;   % maximun bad times per epoch 0.30
P.DefBEa.limBCa       = 0.30;   % maximun channles with some bad data
P.DefBEa.limCCTa      = 1.00;   % maximun interpolated data

P.DefBEdT.limDist     = 2;
P.DefBEdT.limBadDist  = 0.30; %0.30 | 0.10/diff(Epoch.tw);
P.DefBEdT.maxloops    = 1;
P.DefBEdT.rmvmean     = 0 ;

% check the input
[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('eega_avgdata: Non recognized inputs')
end

[Ne, Ns, Nt] = size(EEG.data);
nodata = 0;



%% ------------------------------------------------------------------------
%% Reshape
if ~nodata
    if isempty(P.reshape.nWords) || isempty(P.reshape.tlims) || (P.reshape.nWords==1 && P.reshape.tlims(1)<=EEG.times(1) && P.reshape.tlims(2)>=EEG.times(2))
        P.reshape.do = 0;
    end
    nsmpl = P.reshape.nWords*diff(P.reshape.tlims)/1000*EEG.srate;
    fprintf('### Reshaping the data ###\n')
    [EEG, trlsid] = dft_reshapenooverlap(EEG, P.reshape.nWords, nsmpl, P.reshape.f, P.reshape.fkeep, P.reshape.tlims);
    fprintf('\n')
    [Ne, Ns, Nt] = size(EEG.data);
        
end


%% ------------------------------------------------------------------------
%% Compute BTBC
if ~nodata
    EEG = eega_tDefBTBC(EEG, P.BTep.nbc,P.BCep.nbt,P.BCall.nbt,'keeppre','BT','minBadTime',P.BTep.minBadTime,'minGoodTime',P.BTep.minGoodTime,'maskTime',P.BTep.maskTime);
end

% %% ------------------------------------------------------------------------
% %% Replace the artifacts by...
% if ~nodata
%     EEG = eega_rmvbaddata(EEG, 'BadData', P.artrmv,'TableCNDs',factorCND,'what','BTBC');
% end

%% ------------------------------------------------------------------------
%% Define and remove extra bad epochs
if ~nodata
    P.DefBEa.limBTa = 0.3;
    EEG = eega_tDefBEbaddata(EEG, P.DefBEa,'keeppre',0,'plot', 0);
    EEG = eega_tDefBEdist(EEG, P.DefBEdT.limDist, P.DefBEdT.limBadDist,'maxloops', P.DefBEdT.maxloops,'keeppre',1,'rmvmean',P.DefBEdT.rmvmean, 'plot', 0, 'savefigure', 0);
    EEG = eega_rmvbadepochs(EEG);
    if isempty(EEG.data)
        nodata = 1;
    end
end

%% ------------------------------------------------------------------------
%% Reference average
if ~nodata
    EEG = eega_refavg(EEG);
end

%% ------------------------------------------------------------------------
%% Normalize per long epoch
if ~nodata
    EEG = eega_normalization(EEG,'epochs','single','electrodes','all','mean',0,'sd',[],'BadData','replacebynan');
end


end