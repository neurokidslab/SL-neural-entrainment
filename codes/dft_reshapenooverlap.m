function [EEGout, trialsidx] = dft_reshapenooverlap(EEG, nepochs, samples, factorcond, factors, tlims)

if nargin<6
    tlims = [];
end
mingoodepochs = 0.70;

if ~isempty(EEG.data)
    
    fprintf('### Reshaping the data ###\n')

    %% Take the time
    if ~isempty(tlims)
        %     idxtime = (EEG.times>=tlims(1)) & (EEG.times<=tlims(2));
        idxtime = (EEG.times>=tlims(1)) & (EEG.times<tlims(2));
        EEG.data = EEG.data(:,idxtime,:);
        EEG.times = EEG.times(idxtime);
        if isfield(EEG, 'artifacts')
            EEG.artifacts.BCT = EEG.artifacts.BCT(:,idxtime,:);
            EEG.artifacts.BT = EEG.artifacts.BT(:,idxtime,:);
            EEG.artifacts.CCT = EEG.artifacts.CCT(:,idxtime,:);
        end
    end

    [nCh, nSmp, nEp] = size(EEG.data);

    %% get the condition of each trial
    if ~isempty(factorcond)
        fRESH = cnd_buildtablecond(factorcond, EEG.F);
        [cndsID, cndsVAL, trialsxCND, ~] = cnd_findtrialsxcond(fRESH, EEG.F);
    else
        cndsVAL = {'cnd'};
        cndsID = ones(1,size(EEG.data,3));
    end
    % id_fcnd = false(1,length(EEG.F));
    % for i=1:length(EEG.F)
    %     id_fcnd(i) = strcmp(EEG.F{i}.name,factorcond);
    % end
    % cndsID = EEG.F{id_fcnd}.g;
    % cndsVAL = EEG.F{id_fcnd}.val;

    %% factors
    F = cell(size(factors));
    Fi = cell(size(factors));
    Fg = cell(size(factors));
    Fv = cell(size(factors));
    for ifac=1:length(factors)
        Fi{ifac} = false(1,length(EEG.F));
        for i=1:length(EEG.F)
            Fi{ifac}(i) = strcmp(EEG.F{i}.name,factors{ifac});
        end
        Fg{ifac} = EEG.F{Fi{ifac}}.g;
        Fv{ifac} = EEG.F{Fi{ifac}}.val;
        F{ifac}.name =  EEG.F{Fi{ifac}}.name;
        F{ifac}.val =  Fv{ifac};
        F{ifac}.g =  [];
    end

    D = [];
    BCT = [];
    CCT = [];
    BC = [];
    BT = [];
    BE = [];
    BEm = [];
    trialsidx = nan(nEp,1);
    counter = 0;
    for icnd=1:length(cndsVAL)
        epID =  find(cndsID==icnd);
        n = length(epID);
        NEp_ = floor(n/nepochs);  % number of new epochs in this conditon
        NSmp = nSmp*nepochs;

        if ~isempty(samples) && (NSmp~=samples)
            if NSmp>samples
                n2rmv = NSmp-samples;
                sss = (nSmp:nSmp:nSmp*nepochs);
                if length(sss)<n2rmv
                    error('no ready for this')
                end
            else
                error('no ready for this')
            end
        else
            sss = [];
        end

        % remove data if from the last new epoch there is less than x% of the
        % data, otherwise padd withzeros
        ep2rmv = n - NEp_*nepochs;
        if ep2rmv>(nepochs*mingoodepochs)  % padd with zeros if I need to remove more than nepochs*x epochs
            epoch2add = nepochs-ep2rmv;
            NEp = NEp_+1;
            datapadd = zeros(nCh, nSmp, epoch2add);
            bctpadd = true(nCh, nSmp, epoch2add);
            btpadd = true(1, nSmp, epoch2add);
            bcpadd = false(nCh, 1, epoch2add);
            bepadd = false(1, 1, epoch2add);
        else  % remove extra epochs
            if mod(ep2rmv,2)
                epID = epID(floor(ep2rmv/2)+2:end-floor(ep2rmv/2));
            else
                epID = epID(ep2rmv/2+1:end-ep2rmv/2);
            end
            epoch2add = 0;
            NEp = NEp_;
            datapadd = [];
            bctpadd = [];
            btpadd = [];
            bcpadd = [];
            bepadd = [];
        end

        % take the data and reshape
        if ~isempty(epID)
            % save the position of the original epochs in the new epochs
            id = repmat((counter+1:counter+NEp),[(length(epID)+epoch2add)/NEp 1]);
            id = id(:);
            if epoch2add~=0; id = id(1:end-epoch2add); end
            trialsidx(epID) = id(:);
            counter = counter+NEp;
            % add the padding to the data and reshape the data
            d = cat(3,EEG.data(:,:,epID),datapadd);
            d = reshape(d,[nCh NSmp NEp]);
            % add the padding to the artifacts and reshape the artifacts
            if isfield(EEG, 'artifacts')
                bct = cat(3, EEG.artifacts.BCT(:,:,epID), bctpadd);
                bct = reshape(bct,[nCh NSmp NEp]);
                cct = cat(3, EEG.artifacts.CCT(:,:,epID), bctpadd);
                cct = reshape(cct,[nCh NSmp NEp]);
                bc = cat(3, EEG.artifacts.BC(:,:,epID), bcpadd);
                bc = reshape(bc,[nCh nepochs NEp]);
                bt = cat(3, EEG.artifacts.BT(:,:,epID), btpadd);
                bt = reshape(bt,[1 NSmp NEp]);
                be = cat(3, EEG.artifacts.BE(:,:,epID), bepadd);
                be = reshape(be,[1 nepochs NEp]);
                if isfield(EEG.artifacts, 'BEm')
                    bem = cat(3, EEG.artifacts.BEm(:,:,epID), bepadd);
                    bem = reshape(bem,[1 nepochs NEp]);
                end
            end

            % remove the extra samples
            if ~isempty(sss)
                error('extra samples')
                sss_i = randperm(length(sss));
                sss_i = sss_i(1:n2rmv);

                d(:,sss_i,:) = [];
                if isfield(EEG, 'artifacts')
                    bct(:,sss_i,:) = [];
                    cct(:,sss_i,:) = [];
                    bc(:,sss_i,:) = [];
                    bt(:,sss_i,:) = [];
                    be(:,sss_i,:) = [];
                    if isfield(EEG.artifacts, 'BEm')
                        bem(:,sss_i,:) = [];
                    end
                end
            end

            % Concatenate
            D = cat(3,D,d);
            if isfield(EEG, 'artifacts')
                BCT = cat(3,BCT,bct);
                CCT = cat(3,CCT,cct);
                BC = cat(3,BC,any(bc,2));
                BT = cat(3,BT,bt);
                BE = cat(3,BE,any(be,2));
                if isfield(EEG.artifacts, 'BEm')
                    BEm = cat(3,BEm,any(bem,2));
                end
            end

            % factor
            for ifac=1:length(factors)
                if epoch2add>0
                    gpadd = repmat(Fg{ifac}(epID(end)),[1 epoch2add]);
                    g = cat(2,Fg{ifac}(epID),gpadd);
                    g = reshape(g,[nepochs NEp]);
                else
                    g = reshape(Fg{ifac}(epID),[nepochs NEp]);
                end
                G = nan(1,NEp);
                for ig=1:NEp
                    v = unique(g(:,ig));
                    nv = nan(1,length(v));
                    for iv=1:length(v)
                        nv(iv) = sum(v==iv);
                    end
                    [~, didx] = max(nv);
                    G(ig) = v(didx);
                end
                F{ifac}.g = cat(2,F{ifac}.g,G);
            end
        end

    end


    if isfield(EEG, 'filename')
        [~,EEGout.filename,~] = fileparts(EEG.filename);
    end
    if isfield(EEG, 'filepath')
        EEGout.filepath = EEG.filepath;
    end

    EEGout.data = D;
    if isfield(EEG, 'artifacts')
        EEGout.artifacts.BCT = logical(BCT);
        EEGout.artifacts.CCT = logical(CCT);
        EEGout.artifacts.BC = logical(BC);
        EEGout.artifacts.BT = logical(BT);
        EEGout.artifacts.BE = logical(BE);
        if isfield(EEG.artifacts, 'BEm')
            EEGout.artifacts.BEm = logical(BEm);
        end
        EEGout.artifacts.algorithm.parameters = [];
        EEGout.artifacts.algorithm.stepname = [];
        EEGout.artifacts.algorithm.rejxstep = [];

    end

    EEGout.srate = EEG.srate;
    [EEGout.nbchan, EEGout.pnts, EEGout.trials] = size(D);
    EEGout.times = (1:EEGout.pnts) * 1/EEGout.srate * 1000;
    EEGout.xmin = EEG.times(1);
    EEGout.xmax = EEG.times(end);
    if isfield(EEG, 'chanlocs')
        EEGout.chanlocs = EEG.chanlocs;
        EEGout.chaninfo = EEG.chaninfo;
    end
    EEGout.F = F;

    fprintf('\n')

end
end