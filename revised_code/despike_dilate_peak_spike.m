function [xnew, spkindx, spkinfo] = despike_dilate_peak_spike( ...
    x, spk_thresh, exclusion_range, do_log, spk_side, permission_range, fs)

% Identifies spikes from time-series data;
% removes data and replaces using cubic interpolation.
%
% MINIMAL additions (spike-specific):
%   - Detect on high-passed (or diff) signal to avoid slow artifacts
%   - Optional derivative gate: must have large |dx| too
%   - Reject segments wider than max_width_ms
%
% Added (does NOT affect detection/despiking):
%   - If nargout>=3: compute spike-triggered average (STA) in raw domain
%     and baseline-corrected STA (per epoch, pre-spike baseline window).
%
% INPUT
%   x  = [Nsamples x Nsignals] matrix
%   spk_thresh      = threshold (IQR multiples) for amplitude (default 2.5)
%   exclusion_range = # samples before/after to remove (default 0)
%   do_log          = log-transform before detection (default false)
%   spk_side        = 0 both, -1 negative only, +1 positive only (default 0)
%   permission_range= logical/0-1 vector length Nsamples (optional)
%   fs              = sampling rate in Hz (recommended; used for HP + width + STA)
%
% OUTPUT
%   xnew    = despiked data (raw-domain interpolation)
%   spkindx = linear indices of PEAK sample of each detected spike
%   spkinfo = struct with STA outputs (only if requested)
%
% To convert to time later (single channel case):
%   [samp, chan] = ind2sub(size(x), spkindx);
%   spike_times_s = samp / fs;

if nargin < 7; fs = []; end
if nargin < 6 || isempty(permission_range); permission_range = ones(size(x,1),1); end
if nargin < 5 || isempty(spk_side); spk_side = 0; end
if nargin < 4 || isempty(do_log); do_log = false; end
if nargin < 3 || isempty(exclusion_range); exclusion_range = 0; end
if nargin < 2 || isempty(spk_thresh); spk_thresh = 2.5; end

% ---- defaults (conservative) ----
use_highpass_for_detection = true;
hp_cutoff_hz = 50;
use_derivative_gate = true;
d_thresh = 5;
max_width_ms = 10;

% ---- STA settings (only used if you request spkinfo) ----
sta_win_s = 0.10;                 % +/- around spike
baseline_win_s = [-0.10, -0.02];  % baseline window relative to spike (seconds)

[Nsamp, Nsig] = size(x);
xnew = x;

permission_range = permission_range(:) ~= 0;

% ---------- 0) choose detection signal (keep interpolation on raw x) ----------
xd = x;

if do_log
    xd = log(max(xd, eps));
end

if use_highpass_for_detection
    if isempty(fs)
        xd = [zeros(1,Nsig); diff(xd,1,1)];
    else
        [b,a] = butter(2, hp_cutoff_hz/(fs/2), 'high');
        xd = filtfilt(b,a, xd);
    end
end

% ---------- 1) detect outliers (spike mask) ----------
meds = median(xd, 1, 'omitnan');
iqs  = iqr(xd, 1);
iqs  = max(iqs, eps);

if spk_side
    amp_mask = bsxfun(@gt, spk_side*bsxfun(@minus, xd, meds), spk_thresh * iqs);
else
    amp_mask = bsxfun(@gt, abs(bsxfun(@minus, xd, meds)), spk_thresh * iqs);
end

% ---- derivative gate ----
if use_derivative_gate
    dxd  = [zeros(1,Nsig); diff(xd,1,1)];
    dmed = median(dxd, 1, 'omitnan');
    diq  = iqr(dxd, 1);
    diq  = max(diq, eps);

    d_mask = bsxfun(@gt, abs(bsxfun(@minus, dxd, dmed)), d_thresh * diq);
    spk = amp_mask & d_mask;
else
    spk = amp_mask;
end

spk = bsxfun(@and, spk, permission_range);

% keep original spike mask for peak-finding
spk_orig = spk;

% ---------- 2) interpolate over (dilated) spike windows in RAW x ----------
if any(spk(:))
    for isig = 1:Nsig
        s = spk(:,isig);
        if ~any(s), continue; end

        if exclusion_range > 0
            s = dilate_binary_vector_conv(s, exclusion_range);
        end

        tgood = find(~s);
        tbad  = find(s);

        if numel(tgood) < 2 || isempty(tbad), continue; end

        xnew(tbad,isig) = interp1(tgood, x(tgood,isig), tbad, 'pchip', 'extrap');
    end
end

% ---------- 3) collapse each spike window to a single peak sample ----------
peak_mask = false(Nsamp, Nsig);

for isig = 1:Nsig
    s = spk_orig(:,isig);
    if ~any(s), continue; end

    starts = find(diff([0; s]) == 1);
    ends   = find(diff([s; 0]) == -1);

    sigd = xd(:,isig);

    for k = 1:numel(starts)
        seg_idx = starts(k):ends(k);

        % width gate (avoid slow stuff)
        if ~isempty(fs) && ~isempty(max_width_ms)
            max_w = round((max_width_ms/1000) * fs);
            if numel(seg_idx) > max_w
                continue
            end
        end

        [~, rel_pk] = max(abs(sigd(seg_idx) - meds(isig)));
        pk_samp = seg_idx(rel_pk);

        peak_mask(pk_samp, isig) = true;
    end
end

spkindx = find(peak_mask);


% --- refractory collapse: keep at most one spike per window per channel ---
if ~isempty(fs)
    refrac_ms   = 50;                        % try 50 ms first
    refrac_samp = round(refrac_ms/1000 * fs);

    [spk_samp, spk_chan] = ind2sub([Nsamp, Nsig], spkindx);

    keep = true(size(spk_samp));

    for c = 1:Nsig
        idx = find(spk_chan == c);
        if numel(idx) < 2, continue; end

        [s_sorted, ord] = sort(spk_samp(idx));
        idx = idx(ord);

        last_kept = s_sorted(1);
        for j = 2:numel(idx)
            s0 = spk_samp(idx(j));
            if s0 - last_kept < refrac_samp
                keep(idx(j)) = false;
            else
                last_kept = s0;
            end
        end
    end

    spkindx = spkindx(keep);
end


% ---------- 4) OPTIONAL: compute STA + baseline-corrected STA ----------
spkinfo = struct();
if nargout >= 3
    spkinfo.fs = fs;

    if isempty(fs)
        warning('spkinfo requested but fs is empty. STA time axis will be in samples, not seconds.');
    end

    if isempty(spkindx)
        spkinfo.sta_time = [];
        spkinfo.sta_raw  = [];
        spkinfo.sta_bc   = [];
        spkinfo.note     = 'No spikes detected.';
        return
    end

    % Convert linear peak indices to (sample, channel)
    [spk_samp, spk_chan] = ind2sub([Nsamp, Nsig], spkindx);

    if isempty(fs)
        winS  = round(sta_win_s * 1); % meaningless, but keep safe
        winS  = 0;
        baseS = [0 0];
        t = (-winS:winS).';
    else
        winS  = round(sta_win_s * fs);
        baseS = round(baseline_win_s * fs); % e.g. [-0.1 -0.02] * fs
        t = (-winS:winS).' / fs;
    end

    % Collect epochs from RAW x (not xd)
    L = 2*winS + 1;
    epochs = nan(L, numel(spkindx));

    valid = true(numel(spkindx),1);
    for i = 1:numel(spkindx)
        s0 = spk_samp(i);
        a  = s0 - winS;
        b  = s0 + winS;
        if a < 1 || b > Nsamp
            valid(i) = false;
            continue
        end
        epochs(:,i) = x(a:b, spk_chan(i));
    end

    epochs = epochs(:,valid);
    spkinfo.n_epochs = size(epochs,2);

    if isempty(epochs)
        spkinfo.sta_time = t;
        spkinfo.sta_raw  = [];
        spkinfo.sta_bc   = [];
        spkinfo.note     = 'All epochs were clipped by edges.';
        return
    end

    % Raw STA
    spkinfo.sta_time = t;
    spkinfo.sta_raw  = mean(epochs, 2, 'omitnan');

    % Baseline-correct each epoch using pre-spike baseline window
    if ~isempty(fs)
        base_a = winS + 1 + baseS(1);
        base_b = winS + 1 + baseS(2);
        base_a = max(1, base_a);
        base_b = min(L, base_b);

        if base_a >= base_b
            warning('Baseline window collapsed; skipping baseline correction.');
            spkinfo.sta_bc = spkinfo.sta_raw;
        else
            base = mean(epochs(base_a:base_b,:), 1, 'omitnan');
            epochs_bc = epochs - base; % implicit expansion
            spkinfo.sta_bc = mean(epochs_bc, 2, 'omitnan');
        end
        spkinfo.baseline_win_s = baseline_win_s;
    else
        spkinfo.sta_bc = spkinfo.sta_raw; % no fs -> no baseline window in seconds
        spkinfo.baseline_win_s = [];
    end
end

end


function w = dilate_binary_vector_conv(v, N)
v = v(:) ~= 0;
if N <= 0, w = v; return; end
w = conv(double(v), ones(2*N+1,1), 'same') > 0;
end
