function EEGdft = dft_applyent(EEG, factorCND, P, exclude_bc)

if nargin<4
    exclude_bc = 0;
end

%% ------------------------------------------------------------------------
%% Find the condition of each trial
if ~istable(factorCND)
    tCND = cnd_buildtablecond(factorCND, EEG.F);
else
    tCND = factorCND;
end
[condition, CNDnames, trialsxCND, F] = cnd_findtrialsxcond(tCND, EEG.F);
theCND = unique(condition);
Favg{1}.name = 'CND';
Favg{1}.val = CNDnames';
Favg{1}.g = (1:length(theCND));


% %% ------------------------------------------------------------------------
% %% DSS for denoising
% [Ne, Ns, Nt] = size(EEG.data);
% if ~isempty(DSS)
%     [EEG, ~, ~] = dss_denoise_EEG(EEG, DSS.kpca1, DSS.kpca2, DSS.fbias, DSS.fapply, 'timewind', DSS.timewind );
% end

%% ------------------------------------------------------------------------
%% Compute the entreinment
[Ne, Ns, Nt] = size(EEG.data);
if ~isempty(EEG.data)
    
    %% FFT
    fprintf('Computing FFT and entreinment...\n')
    Y = EEG.data;
    Y = Y - repmat(mean(Y,2), [1 Ns 1]); 
    f_res = EEG.srate / Ns;             % resolution (Hz)
    freq  = (0:floor(Ns/2))*f_res;  % frequencies
    Nfrq = length(freq);
    Ydft = fft(Y,Ns,2);
    Ydft = Ydft(:,1:floor(Ns/2)+1,:);
    W = repmat(EEG.artifacts.BC==0, [1 size(Ydft,2) 1]); % weights to reject bad electrodes per epoch
    
    %% PLV per condition
    Ncnd = length(unique(condition));
    data_plv = nan(Ne,Nfrq,Ncnd);
    for i=1:Ncnd
        if exclude_bc==0
            data_plv(:,:,i) = abs( mean( Ydft(:,:,condition==i) ./ abs(Ydft(:,:,condition==i)) , 3 ) );
        else
            m = Ydft(:,:,condition==i) ./ abs(Ydft(:,:,condition==i));
            n = sum(W(:,:,condition==i),3);
            data_plv(:,:,i) = abs( sum( m .* W(:,:,condition==i)  , 3 )./n );
        end
    end
    
    % SNR
    if isfield(P,'typeSNRPWR') && ~strcmp(P.typeSNRPLV, 'none')
        fprintf('- PLV: SNR %s...\n',P.typeSNRPLV)
        data_plv = dft_norm( data_plv, freq, 'binsnorm', P.binsSNR,'fband',P.fband,'ftarget',P.ftarget,...
            'typeSNR',P.typeSNRPLV,'noutliers',P.noutSNR,'silent',1);
    end
    
    % Z-score to be sure it is centered at zero
    if isfield(P,'zscorePWR') && P.zscorePLV
        fprintf('- PLV: Z-score...\n')
        data_plv = dft_norm( data_plv, freq, 'fband', P.foutput, 'ftarget', P.ftarget,...
            'binsnorm','all','typeSNR','zscore','noutliers',P.noutSNR,'silent',1);
    end
    
    %% Power Spectrum per condition
    Ncnd = length(unique(condition));
    data_pwr = nan(Ne,Nfrq,Ncnd);
    for i=1:Ncnd
        if exclude_bc==0
            ydftavg = mean( Ydft(:,:,condition==i), 3);
            data_pwr(:,:,i) = single( ydftavg .* conj(ydftavg) );  % Power Spectrum
        else
            n = sum(W(:,:,condition==i),3);
            ydftavg = sum( Ydft(:,:,condition==i) .* W(:,:,condition==i), 3)./n;
            data_pwr(:,:,i) = single( ydftavg .* conj(ydftavg) );  % Power Spectrum
        end
    end

    % SNR
    if isfield(P,'typeSNRPWR') && ~strcmp(P.typeSNRPWR, 'none')
        fprintf('- PWR: SNR %s...\n',P.typeSNRPWR)
        data_pwr = dft_norm( data_pwr, freq, 'binsnorm', P.binsSNR,'fband',P.fband,'ftarget',P.ftarget,...
        'typeSNR',P.typeSNRPWR,'noutliers',P.noutSNR,'silent',1);
    end
    
    % Z-score to be sure it is centered at zero
    if isfield(P,'zscorePWR') && P.zscorePWR
        fprintf('- PWR: Z-score...\n')
        data_pwr = dft_norm( data_pwr, freq, 'fband', P.foutput, 'ftarget', P.ftarget,...
            'binsnorm','all','typeSNR','zscore','noutliers',P.noutSNR,'silent',1);
    end
    
    
    %% Keep the freqeuncies in the output frequency range
    if isfield(P,'foutput') && ~isempty(P.foutput)
        idx_freq_out = freq>=P.foutput(1) & freq<=P.foutput(end);
        freq = freq(idx_freq_out);
        data_plv = data_plv(:,idx_freq_out,:);
        data_pwr = data_pwr(:,idx_freq_out,:);
    end
    
else
    data_plv = [];
    data_pwr = [];
    freq = [];
end

%% ------------------------------------------------------------------------
%% Store the data
EEGdft.filename         = EEG.filename;
EEGdft.TrialsxCND       = array2table(trialsxCND,'VariableNames',CNDnames');
EEGdft.chaninfo         = EEG.chaninfo;
EEGdft.chanlocs         = EEG.chanlocs;
EEGdft.srate            = EEG.srate;
EEGdft.times            = EEG.times;
EEGdft.freqs            = freq;
EEGdft.F                = Favg;
if isfield(EEG,'covnoise')
    EEGdft.covnoise = EEG.covnoise;
end
EEGdft.plv             = data_plv;
EEGdft.pwr             = data_pwr;
    

end