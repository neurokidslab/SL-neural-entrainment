function EEGout = dft_reshape(EEG, nepochs, factorcond, factors, tlims)

if ~isempty(EEG.data)

    nsmpl = nepochs*diff(tlims)/1000*EEG.srate;

    [EEGout, trialsidx] = dft_reshapenooverlap(EEG, nepochs, nsmpl, factorcond, factors, tlims);
else 
    EEGout = EEG;
end
end
