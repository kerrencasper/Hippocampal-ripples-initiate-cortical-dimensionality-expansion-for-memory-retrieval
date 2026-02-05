function alldat = remove_ripple_duplicates(alldat)



all_ripples = [];

for ichannel = 1:numel(alldat)
    all_ripples = [all_ripples alldat{ichannel}.evtIndiv.maxTime];
end

unique_ripples = unique(all_ripples);

ripples_cooccur = [];

for ii=1:numel(unique_ripples)
    if sum(all_ripples==unique_ripples(ii)) > 1
        ripples_cooccur = [ripples_cooccur unique_ripples(ii)];
    end
end


for ichannel = 1:size(alldat,1)

    evs = alldat{ichannel}.evtIndiv.maxTime;

    these_ripples_duplicates = ismember(evs,ripples_cooccur);


    alldat{ichannel}.evtIndiv.staTime(these_ripples_duplicates)      = [];
    alldat{ichannel}.evtIndiv.midTime(these_ripples_duplicates)      = [];
    alldat{ichannel}.evtIndiv.endTime(these_ripples_duplicates)      = [];
    alldat{ichannel}.evtIndiv.stage(these_ripples_duplicates)        = [];
    alldat{ichannel}.evtIndiv.duration(these_ripples_duplicates)     = [];
    alldat{ichannel}.evtIndiv.maxTime(these_ripples_duplicates)      = [];
    alldat{ichannel}.evtIndiv.minTime(these_ripples_duplicates)      = [];
    alldat{ichannel}.evtIndiv.minAmp(these_ripples_duplicates)       = [];
    alldat{ichannel}.evtIndiv.maxAmp(these_ripples_duplicates)       = [];
    alldat{ichannel}.evtIndiv.envMaxAmp(these_ripples_duplicates)    = [];
    alldat{ichannel}.evtIndiv.envMaxTime(these_ripples_duplicates)   = [];
    alldat{ichannel}.evtIndiv.envMean(these_ripples_duplicates)      = [];
    alldat{ichannel}.evtIndiv.envSum(these_ripples_duplicates)       = [];
    alldat{ichannel}.evtIndiv.peaks(these_ripples_duplicates)        = [];
    alldat{ichannel}.evtIndiv.troughs(these_ripples_duplicates)      = [];
    alldat{ichannel}.evtIndiv.freq(these_ripples_duplicates)         = [];

    alldat{ichannel}.evtIndiv.numEvt = sum(not(these_ripples_duplicates));

end

these_ripples_duplicates = [];


end