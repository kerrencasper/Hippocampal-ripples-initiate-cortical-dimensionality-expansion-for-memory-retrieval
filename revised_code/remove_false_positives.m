function alldat = remove_false_positives(alldat)


for ichannel=1:numel(alldat)

    false_positives     = alldat{ichannel}.falseposRjct.rejects{1};

    alldat{ichannel}.evtIndiv.staTime(false_positives)      = [];
    alldat{ichannel}.evtIndiv.midTime(false_positives)      = [];
    alldat{ichannel}.evtIndiv.endTime(false_positives)      = [];
    alldat{ichannel}.evtIndiv.stage(false_positives)        = [];
    alldat{ichannel}.evtIndiv.duration(false_positives)     = [];
    alldat{ichannel}.evtIndiv.maxTime(false_positives)      = [];
    alldat{ichannel}.evtIndiv.minTime(false_positives)      = [];
    alldat{ichannel}.evtIndiv.minAmp(false_positives)       = [];
    alldat{ichannel}.evtIndiv.maxAmp(false_positives)       = [];
    alldat{ichannel}.evtIndiv.envMaxAmp(false_positives)    = [];
    alldat{ichannel}.evtIndiv.envMaxTime(false_positives)   = [];
    alldat{ichannel}.evtIndiv.envMean(false_positives)      = [];
    alldat{ichannel}.evtIndiv.envSum(false_positives)       = [];
    alldat{ichannel}.evtIndiv.peaks(false_positives)        = [];
    alldat{ichannel}.evtIndiv.troughs(false_positives)      = [];
    alldat{ichannel}.evtIndiv.freq(false_positives)         = [];


    alldat{ichannel}.evtIndiv.numEvt = sum(not(false_positives));
end

end