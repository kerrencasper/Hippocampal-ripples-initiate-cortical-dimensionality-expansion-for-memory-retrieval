function [inData,alldat] = detect_ripples(data, settings)

alldat = cell(numel(data.label),1);

for ichannel = 1:numel(data.label)

    inData              = {};
    inData.label        = data.label(ichannel);
    inData.trial        = {data.trial{1}(ichannel,:)};
    inData.time         = data.time;
    inData.fsample      = data.fsample;
    inData.sampleinfo   = [1 numel(data.time{1})];

    inData.artifacts    = {data.artifact{ichannel}};
    inData.staging      = {zeros(size(data.time{1}))};


    cfg                 = [];
    cfg.bsfilter        = 'yes';
    cfg.bsfreq          = [49 51;99 101;149 151;199 201];
    cfg.detrend         = 'yes';
    tmp                 = ft_preprocessing(cfg,inData);

    inData.trial        = tmp.trial;
    tmp = [];

    tfg             = [];
    tfg.thresCh     = {'all'};
    tfg.detectCh    = {'all'};

    tfg.param.bpfreq    = [80 120]; %80 120
    tfg.param.stageoi   = 0;
    tfg.param.artfctPad = [-0.05 0.05]; % plus minus 50 ms.
    tfg.param.thresType = 'channelwise';
    tfg.param.filtType  = 'fir';
    tfg.param.envType   = 'rms';
    tfg.param.envWin    = 0.02;

    tfg.criterion.len       = [0.038 0.5];
    tfg.criterion.center    = 'mean';
    tfg.criterion.var       = 'centerSD';
    tfg.criterion.scale     = 1.5; %2.5

    tfg.paramOpt.smoothEnv   = 0.02;
    tfg.paramOpt.upperCutoff = 9;

    tfg.paramOpt.minCycles       = 1;
    tfg.paramOpt.minCyclesNum    = 3;

    tfg.doFalseposRjct          = settings.remove_falsepositives;
    tfg.falseposRjct.freqlim    = [75 125];
    tfg.falseposRjct.timePad    = [-0.25 0.25];
    tfg.falseposRjct.tfrFreq    = 65 : 2 : 135;
    tfg.falseposRjct.tfrTime    = -0.1 : 0.002 : 0.1;
    tfg.falseposRjct.tfrWin     = ceil(0.1 * tfg.falseposRjct.tfrFreq) ./ tfg.falseposRjct.tfrFreq;
    tfg.falseposRjct.avgWin     = [-0.05 0.05];
    tfg.falseposRjct.prominence = 0.2;

    outRipple = Slythm_DetectSpindles_v2(tfg, inData);

    alldat{ichannel} = outRipple;

end

end