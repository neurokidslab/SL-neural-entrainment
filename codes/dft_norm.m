% INPUTS
% EEG       : data
% fband     :

% OUTPUTS
% EEGf    :

function [ PSn ] = dft_norm( PS , freq, varargin )

%% ------------------------------------------------------------------------
%% Parameters

% Default parameters
P.fband = [];
P.ftarget = [];
P.noutliers = 2.0;  % number of standar deviations used to define outliers, which won't be used for normalization
P.binsnorm = [(-11:-2) (2:11)];% 10 bins before and 10 after but not the inmidiatlly attached
P.typeSNR = 'zscore';
P.silent = 0;

[P, OK, extrainput] = eega_getoptions(P, varargin);
if ~OK
    error('dft_powernorm: Non recognized inputs')
end

if isempty(P.fband)
    P.fband =[freq(1) freq(end)];
end

if ischar(P.binsnorm) && strcmp(P.binsnorm,'all')
    doall = 1;
else
    doall = 0;
end
%% ------------------------------------------------------------------------
%% Normalize the data

if ~P.silent
    fprintf(sprintf('Normalization\n'))
end

% find the indexes of the target frequencies to exclude
idxftarget=nan(1,length(P.ftarget));
for i=1:length(P.ftarget)
    [~,idxftarget(i)]=min(abs(freq-P.ftarget(i)));
end
ftarget=freq(idxftarget);

% normalize
PSn = nan(size(PS));
if doall
    % find the indexes of the frequencies for normalization
    f_idx = (freq>=P.fband(1)) & (freq<=P.fband(end)) & ~ismember(freq, ftarget);

    switch P.typeSNR
        case 'zscore'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    f_idx_f = f_idx;
                    y = PS(el,f_idx_f,i);
                    y = y(:);
                    y = removeoutliers(y, P.noutliers);
                    mu = nanmean(y);
                    sd = nanstd(y);
                    PSn(el,:,i) = (PS(el,:,i) - mu) / sd;
                end
            end
        case 'mean'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    f_idx_f = f_idx;
                    y = PS(el,f_idx_f,i);
                    y = y(:);
                    y = removeoutliers(y, P.noutliers);
                    mu = nanmean(y);
                    PSn(el,:,i) = (PS(el,:,i) - mu);
                end
            end
        case 'ratio'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    f_idx_f = f_idx;
                    y = PS(el,f_idx_f,i);
                    y = y(:);
                    y = removeoutliers(y, P.noutliers);
                    mu = nanmean(y);
                    PSn(el,:,i) = (PS(el,:,i) / mu);
                end
            end
        case 'diffratio'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    f_idx_f = f_idx;
                    y = PS(el,f_idx_f,i);
                    y = y(:);
                    y = removeoutliers(y, P.noutliers);
                    mu = nanmean(y);
                    PSn(el,:,i) = (PS(el,:,i) - mu) / mu;
                end
            end
        case 'linear'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    f_idx_f = f_idx;
                    x = freq(f_idx_f);
                    x = x(:);
                    y = PS(el,f_idx_f,i);
                    y = y(:);
                    [x, y] = removeoutlierslfit(x, y, P.noutliers);
                    p = polyfit(x,y,1);
                    Y = polyval(p,freq);
                    PSn(el,:,i) = PS(el,:,i) - Y;
                end
            end
        case 'linearzscore'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    f_idx_f = f_idx;
                    x = freq(f_idx_f);
                    x = x(:);
                    y = PS(el,f_idx_f,i);
                    y = y(:);
                    [x, y] = removeoutlierslfit(x, y, P.noutliers);
                    p = polyfit(x,y,1);
                    Y = polyval(p,freq);
                    noise = std(PS(el,:,i) - Y);
                    PSn(el,:,i) = (PS(el,:,i) - Y) ./ noise;
                end
            end
    end
else
    % find the indexes of the frequencies to use for normalization
    f_idx = (freq>=P.fband(1)) & (freq<=P.fband(end)) & ~ismember(freq, ftarget);
    
    % find the indexes of the frequencies to normalizae
    f2norm_idx = (freq>=P.fband(1)) & (freq<=P.fband(end));
    f2norm = find(f2norm_idx);
    
    switch P.typeSNR
        case 'zscore'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    for fi=1:length(f2norm)
                        f = f2norm(fi);
                        binsnorm = findbinsnorm(freq, P.binsnorm, f);
                        f_idx_f = binsnorm & freq~=freq(f) & f_idx;
                        if sum(f_idx_f)>2
                            y = PS(el,f_idx_f,i);
                            y = y(:);
                            y = removeoutliers(y, P.noutliers);
                            mu = nanmean(y);
                            sd = nanstd(y);
                            PSn(el,f,i) = (PS(el,f,i) - mu) / sd;
                        end
                    end
                end
            end
        case 'mean'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    for fi=1:length(f2norm)
                        f = f2norm(fi);
                        binsnorm = findbinsnorm(freq, P.binsnorm, f);
                        f_idx_f = binsnorm & freq~=freq(f) & f_idx;
                        if sum(f_idx_f)>2
                            y = PS(el,f_idx_f,i);
                            y = y(:);
                            y = removeoutliers(y, P.noutliers);
                            mu = nanmean(y);
                            PSn(el,f,i) = (PS(el,f,i) - mu);
                        end
                    end
                end
            end
        case 'ratio'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    for fi=1:length(f2norm)
                        f = f2norm(fi);
                        binsnorm = findbinsnorm(freq, P.binsnorm, f);
                        f_idx_f = binsnorm & freq~=freq(f) & f_idx;
                        if sum(f_idx_f)>2
                            y = PS(el,f_idx_f,i);
                            y = y(:);
                            y = removeoutliers(y, P.noutliers);
                            mu = nanmean(y);
                            PSn(el,f,i) = (PS(el,f,i) / mu);
                        end
                    end
                end
            end
        case 'logratio'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    for fi=1:length(f2norm)
                        f = f2norm(fi);
                        binsnorm = findbinsnorm(freq, P.binsnorm, f);
                        f_idx_f = binsnorm & freq~=freq(f) & f_idx;
                        if sum(f_idx_f)>2
                            y = PS(el,f_idx_f,i);
                            y = y(:);
                            y = removeoutliers(y, P.noutliers);
                            mu = nanmean(y);
                            PSn(el,f,i) = log10(PS(el,f,i) / mu);
                        end
                    end
                end
            end
        case 'diffratio'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    for fi=1:length(f2norm)
                        f = f2norm(fi);
                        binsnorm = findbinsnorm(freq, P.binsnorm, f);
                        f_idx_f = binsnorm & freq~=freq(f) & f_idx;
                        if sum(f_idx_f)>2
                            y = PS(el,f_idx_f,i);
                            y = y(:);
                            y = removeoutliers(y, P.noutliers);
                            mu = nanmean(y);
                            PSn(el,f,i) = (PS(el,f,i) - mu) / mu;
                        end
                    end
                end
            end
        case 'linear'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    for fi=1:length(f2norm)
                        f = f2norm(fi);
                        binsnorm = findbinsnorm(freq, P.binsnorm, f);
                        f_idx_f = binsnorm & freq~=freq(f) & f_idx;
                        if sum(f_idx_f)>2
                            x = freq(f_idx_f);
                            x = x(:);
                            y = PS(el,f_idx_f,i);
                            y = y(:);
                            [x, y] = removeoutlierslfit(x, y, P.noutliers);
                            p = polyfit(x,y,1);
                            Y = polyval(p,freq(f));
                            PSn(el,f,i) = PS(el,f,i) - Y;
                        end
                    end
                end
            end
        case 'linearzscore'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    for fi=1:length(f2norm)
                        f = f2norm(fi);
                        binsnorm = findbinsnorm(freq, P.binsnorm, f);
                        f_idx_f = binsnorm & freq~=freq(f) & f_idx;
                        if sum(f_idx_f)>2
                            x = freq(f_idx_f);
                            x = x(:);
                            y = PS(el,f_idx_f,i);
                            y = y(:);
                            [x, y] = removeoutlierslfit(x, y, P.noutliers);
                            p = polyfit(x,y,1);
                            Yf = polyval(p,freq(f));
                            Y = polyval(p,freq(f_idx_f));
                            noise = std(PS(el,f_idx_f,i) - Y);
                            PSn(el,f,i) = (PS(el,f,i) - Yf) / noise;
                        end
                    end
                end
            end
            
        case 'logzscore'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    for fi=1:length(f2norm)
                        f = f2norm(fi);
                        binsnorm = findbinsnorm(freq, P.binsnorm, f);
                        f_idx_f = binsnorm & freq~=freq(f) & f_idx;
                        if sum(f_idx_f)>2
                            x = freq(f_idx_f);
                            y = PS(el,f_idx_f,i);
                            y = removeoutliers(y, P.noutliers);
                            mu = nanmean(log10(y));
                            sd = nanstd(log10(y));
                            PSn(el,f,i) = (log10( PS(el,f,i) ) - mu) ./ sd;
                        end
                    end
                end
            end
            
        case 'logmean'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    for fi=1:length(f2norm)
                        f = f2norm(fi);
                        binsnorm = findbinsnorm(freq, P.binsnorm, f);
                        f_idx_f = binsnorm & freq~=freq(f) & f_idx;
                        if sum(f_idx_f)>2
                            x = freq(f_idx_f);
                            y = PS(el,f_idx_f,i);
                            y = removeoutliers(y, P.noutliers);
                            PSn(el,f,i) = 20*log10( PS(el,f,i) ./ mean(y,2) );
                        end
                    end
                end
            end
            
        case 'loglinear'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    for fi=1:length(f2norm)
                        f = f2norm(fi);
                        binsnorm = findbinsnorm(freq, P.binsnorm, f);
                        f_idx_f = binsnorm & freq~=freq(f) & f_idx;
                        if sum(f_idx_f)>2
                            x = log10(freq(f_idx_f));
                            y = log10(PS(el,f_idx_f,i));
                            [x, y] = removeoutlierslfit(x, y, P.noutliers);
                            p = polyfit(x,y,1);
                            Y = polyval(p,log10(freq(f)));
                            PSn(el,f,i) = 20 * (log10( PS(el,f,i) ) - Y );
                        end
                    end
                end
            end
            
        case 'loglinearzscore'
            for i=1:size(PS,3)
                for el=1:size(PS,1)
                    for fi=1:length(f2norm)
                        f = f2norm(fi);
                        binsnorm = findbinsnorm(freq, P.binsnorm, f);
                        f_idx_f = binsnorm & freq~=freq(f) & f_idx;
                        if sum(f_idx_f)>2
                            x = log10(freq(f_idx_f));
                            y = log10(PS(el,f_idx_f,i));
                            [x, y] = removeoutlierslfit(x, y, P.noutliers);
                            p = polyfit(x,y,1);
                            Y = polyval(p,log10(freq(f)));
                            Ynoise = polyval(p,log10(freq(f_idx_f)));
                            noise = std( log10(PS(el,f_idx_f,i)) - Ynoise);
                            PSn(el,f,i) = (log10( PS(el,f,i) ) - Y ) / noise;
                        end
                    end
                end
            end
            
            
    end
end

if ~P.silent; fprintf('\n'); end

end

function binsnorm = findbinsnorm(freq, nbins, f)
binsnorm = zeros(size(freq));
idxbinsnorm = nbins + f;
idxbinsnorm(idxbinsnorm<1)=[];
idxbinsnorm(idxbinsnorm>length(freq))=[];
% binsup = idxbinsnorm>f;
% binsdown = idxbinsnorm<f;
% if sum(binsup)~=sum(binsdown)
%     if sum(binsdown)==0
%         if length(idxbinsnorm)>1
%             idxbinsnorm = idxbinsnorm(1:2);
%         end
%     elseif sum(binsup)==0
%         if length(idxbinsnorm)>1
%             idxbinsnorm = idxbinsnorm(end-1:end);
%         end
%     elseif sum(binsup)<sum(binsdown)
%         n = sum(binsup);
%         idxbinsnorm(1:sum(binsdown)-n) = [];
%     elseif sum(binsup)>sum(binsdown)
%         n = sum(binsdown);
%         idxbinsnorm(2*n+1:end) = [];
%     end
% end
binsnorm(idxbinsnorm) = 1;

end

function [y] = removeoutliers(y, n)

mu = mean(y);
sigma = std(y);
thrs = [mu-n*sigma mu+n*sigma];

idx = y<thrs(1)  | y>thrs(2);
y = y(~idx);

end

function [x, y] = removeoutliers2(x, y, n)

mu = mean(y);
sigma = std(y);
thrs = [mu-n*sigma mu+n*sigma];

idx = y<thrs(1)  | y>thrs(2);
y = y(~idx);
x = x(~idx);

end

function [x, y] = removeoutlierslfit(x, y, n)

p = polyfit(x,y,1);
Y = polyval(p,x);
res = Y-y;

mu = mean(res);
sigma = std(res);
thrs = [mu-n*sigma mu+n*sigma];

idx = res<thrs(1)  | res>thrs(2);
y = y(~idx);
x = x(~idx);

end