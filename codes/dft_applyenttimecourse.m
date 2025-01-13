function EEGdft = dft_applyenttimecourse(EEG, nepxtime, nepstep, varargin)

%% ------------------------------------------------------------------------
%% Parameters

% Remove artifacts
P.exclude_bc        = 0;

% Reshape
P.reshape.f         = [];
P.reshape.fkeep     = [];
P.reshape.nWords    = [];
P.reshape.tlims     = [];

% Define Bad Segments and Electrodes
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


P.mintrl        = 2;

% If DSS is applied, then P.DSS has to have the fields
%       - P.DSS.reshape   
%       - P.DSS.kDSS      
%       - P.DSS.nDSS      
%       - P.DSS.fDSSbias  
%       - P.DSS.fDSSapply 
%       - P.DSS.twDSS
P.DSS.do        = 0;

P.ent.foutput       = [0.2 15];
P.ent.fband         = [];
P.ent.ftarget       = [];
P.ent.binsSNR       = [(-5:-1) (1:5)];
P.ent.noutSNR       = 2.0;

P.ent.typeSNRPLV    = 'zscore';   % 'linearglobal' / 'linearlocal' / 'meanlocal'
P.ent.typeSNRPWR    = 'loglinearzscore';   % 'linearglobal' / 'linearlocal' / 'meanlocal'
P.ent.zscorePLV     = 0;
P.ent.zscorePWR     = 0;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('dft_applyplv: Non recognized inputs')
end

[Ne, Ns, Nt] = size(EEG.data);


if Nt>0
    
    
    %% ------------------------------------------------------------------------
    %% Get relevant info regarding the time slots
    % different epochs
    fEPOCH = cnd_buildtablecond({'Epoch'}, EEG.F);
    [EPOCH, theEPOCH, ~, ~] = cnd_findtrialsxcond(fEPOCH, EEG.F);
    fCND = cnd_buildtablecond({'eventCnd'}, EEG.F);
    [CND, theCND, ~, ~] = cnd_findtrialsxcond(fCND, EEG.F);
    thecnd = find(strcmp(theCND,'STR'));
    % epochs to considered at each time
    TOTep = length(EEG.F{end}.val); % total possible epochs
    Iep = (1:nepstep:TOTep-(nepxtime*P.reshape.nWords)+1);
    ntimes = length(Iep);
    % Create the variables to hold the data
    times = cell(1,ntimes);
    PLV = nan(Ne, 200, ntimes);
    Pwr = nan(Ne, 200, ntimes);
    trialsxTIME = zeros(1,ntimes);
    % create the factor
    Ftime.name = 'TIME';
    Ftime.val = cell(1,ntimes);
    Ftime.g = (1:ntimes);
    for itime=1:ntimes
        
        times{itime} = sprintf('t%03d',itime);
        fprintf('%s\nCOMPUTING AT TIME  %s\n\n',repmat('-',[1 30]),times{itime})
        
        EEGi = EEG;
        Ftime.val{itime} = times{itime};
        [Ne, Ns, Nt] = size(EEG.data);
        
        % Recover the epochs for that time slot
        idxep = (Iep(itime):Iep(itime)+(nepxtime*P.reshape.nWords)-1);
        epval = cell(1,length(idxep));
        for k=1:length(idxep)
            epval{k} = sprintf('x%04d',idxep(k));
        end
        idxEPOCH = ismember(theEPOCH(EPOCH),epval);
        
        % Create a factor defining for each epoch if it is in the slot
        Ftsl.name = 'timeslot';
        Ftsl.val = {'in' 'out'};
        Ftsl.g = 2*ones(1,Nt);
        Ftsl.g(idxEPOCH) = 1;
        EEGi.F = cat(1,EEG.F(1:end-1),Ftsl,EEG.F(end));
        
        % Take the corresponding subset of the data
        EEGi.artifacts.BE(~idxEPOCH) = 1;
        EEGi = eega_rmvbadepochs(EEGi,'Silent',1);
        
        % All the stepd before
        P.reshape_i.nWords = P.reshape.nWords;
        P.reshape_i.f = {'timeslot'};
        P.reshape_i.fkeep = {'timeslot' 'eventCnd'};
        P.reshape_i.tlims = P.reshape.tlims; 
        cnds = [];
        cnds{1} = EEGi.F{1}.g==1;
        cnds{2} = EEGi.F{1}.g==2;
        gblck = [1 3 4 5 6 7 8 9 10 2];
        for j=1:10
            cnds{2+j} = EEGi.F{1}.g==3 & EEGi.F{2}.g==gblck(j);
        end
        [EEGi] = dft_2doforfft(EEGi, 'reshape',P.reshape_i,...
            'BCall',P.BCall,'BCep',P.BCep,'BTep',P.BTep,'DefBEa',P.DefBEa,'DefBEdT',P.DefBEdT);
        trialsxTIME(itime) = size(EEGi.data,3);
        
        if size(EEGi.data,3)>=P.mintrl
            
            [Ne, Ns, Nt] = size(EEGi.data);
            
            %% DSS for denoising
            if P.DSS.do
                % Apply it per time window
                thetimeidx = false(length(EEGi.F),1);
                for ifi=1:length(EEGi.F)
                    thetimeidx(ifi) = strcmp('timeslot',EEGi.F{ifi}.name);
                end
                DSSi = P.DSS;
                mm = { {'timeslot'} {'in'}};
                DSSi.fbias = mm;
                DSSi.fapply = mm;
                [EEGi, ~, ~] = dss_denoise_EEG(EEGi, DSSi.kpca1, DSSi.kpca2, DSSi.fbias, DSSi.fapply, 'timewind', P.DSS.timewind );               
            end
            
            %% Entreiment computation
            EEGidft = dft_applyent(EEGi, {'timeslot'}, P.ent, P.exclude_bc);
            Nfrq = size(EEGidft.plv,2);
            freq = EEGidft.freqs; 
            PLV(:,1:Nfrq,itime) = EEGidft.plv;
            Pwr(:,1:Nfrq,itime) = EEGidft.pwr;
                
        
        else
            fprintf('Not enoutgh good trials\n')
        end
        fprintf('\n')
                
    end
    PLV = PLV(:,1:Nfrq,:);
    Pwr = Pwr(:,1:Nfrq,:);
    
     
    %% Keep the freqeuncies in the output frequency range
    if ~isempty(P.ent.foutput)
        idx_freq_out = freq>=P.ent.foutput(1) & freq<=P.ent.foutput(end);
        freq = freq(idx_freq_out);
        PLV = PLV(:,idx_freq_out,:);
        Pwr = Pwr(:,idx_freq_out,:);
    end
    
else
    Ftime = [];
    PLV = [];
    Pwr = [];
    freq = [];
    
end


%% Store the data
x = size(EEG.data,2)/EEG.srate;
nep = size(EEG.data,3);
time = (x*(1:nepstep:nep) + P.reshape.nWords*nepxtime*x/2 -x)';
time = time(1:ntimes);

EEGdft.filename     = EEG.filename;
EEGdft.chaninfo     = EEG.chaninfo;
EEGdft.chanlocs     = EEG.chanlocs;
EEGdft.srate        = EEG.srate;
EEGdft.freqs        = freq;
EEGdft.times        = time;
EEGdft.F            = {Ftime};
EEGdft.plv          = PLV;
EEGdft.pwr          = Pwr;
EEGdft.TrialsxCND   = array2table(trialsxTIME,'VariableNames',times');


end