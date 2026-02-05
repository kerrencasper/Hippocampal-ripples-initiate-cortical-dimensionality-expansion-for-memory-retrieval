function [xnew, spkindx] = despike_dilate_peak_spike(x, spk_thresh, exclusion_range, do_log, spk_side, permission_range)

% Identifies spikes from time-series data;
% removes data and replaces using cubic interpolation.
%
% INPUT
%   x  = [Nsamples x Nsignals] matrix
%   spk_thresh      = threshold (IQR multiples)
%   exclusion_range = # samples before/after to remove
%   do_log          = log-transform before detection
%   spk_side        = 0 both, -1 negative only, +1 positive only
%   permission_range= logical/0-1 vector of length Nsamples (optional)
%
% OUTPUT
%   xnew    = despiked data
%   spkindx = **peak sample indices** of each detected spike
%
% To convert to time later:  spike_times = spkindx / fs;

if nargin < 6; permission_range = ones(size(x,1),1); end
if nargin < 5; spk_side = 0; end
if nargin < 4; do_log = false; end
if nargin < 3; exclusion_range = 0; end
if nargin < 2; spk_thresh = 2.5; end

[Nsamp, Nsig] = size(x);
xnew = x;

% ---------- 1) detect outliers (spike mask) ----------
if do_log
    y    = log(x);
    meds = median(y);
    iqs  = iqr(y);
    if spk_side
        spk = bsxfun(@gt, spk_side*bsxfun(@minus, y, meds), spk_thresh * iqs);
    else
        spk = bsxfun(@gt, abs(bsxfun(@minus, y, meds)), spk_thresh * iqs);
    end
else
    meds = median(x);
    iqs  = iqr(x);
    if spk_side
        spk = bsxfun(@gt, spk_side*bsxfun(@minus, x, meds), spk_thresh * iqs);
    else
        spk = bsxfun(@gt, abs(bsxfun(@minus, x, meds)), spk_thresh * iqs);
    end
end

spk = bsxfun(@and, spk, permission_range);  % apply permission mask

% keep a copy of the original spike mask for peak-finding
spk_orig = spk;

% ---------- 2) interpolate over (dilated) spike windows ----------
if any(spk(:))   % if we find any spikes at all
    for isig = 1:Nsig
        s = spk(:,isig);            % spike mask for this channel
        if any(s)
            if exclusion_range > 0
                % expand spike regions for interpolation
                s = dilate_binary_vector(s, exclusion_range);
            end

            tgood = find(~s);
            tbad  = find(s);

            xnew(tbad,isig) = interp1(tgood, x(tgood,isig), tbad, 'pchip');
        end
    end
end

% ---------- 3) collapse each spike window to a single peak sample ----------
peak_mask = false(Nsamp, Nsig);

for isig = 1:Nsig
    s = spk_orig(:,isig);   % use ORIGINAL mask (without dilation) for peaks
    if ~any(s), continue; end

    % find contiguous segments of spike samples
    starts = find(diff([0; s]) == 1);    % 0 -> 1
    ends   = find(diff([s; 0]) == -1);   % 1 -> 0

    sig = x(:,isig);  % original data for this channel

    for k = 1:numel(starts)
        seg_idx = starts(k):ends(k);

        % peak = sample with max |deviation| within this segment
        [~, rel_pk] = max(abs(sig(seg_idx) - meds(isig)));
        pk_samp     = seg_idx(rel_pk);

        peak_mask(pk_samp, isig) = true;
    end
end

% Return linear indices of **peak samples only**
spkindx = find(peak_mask);

end


function w = dilate_binary_vector(v, N)
    % dilates blocks of ones in a binary vector by N samples
    v = double(v); 
    b = ones(N+1,1);
    w = (filtfilt(b,1,v) > 0);
end

