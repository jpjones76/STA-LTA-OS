function [Xf,X0] = genprep(X, r, fs)
% [Xf,X0] = genprep(X,r,fs);
%
% Generic preprocessing of input data X using ruleset r.
%
% INPUTS
% X     Seismic data in column vectors
% r     Ruleset as a cell string. See below.
% fs    Sampling frequency.
%
% OUTPUTS
% Xf    Processed time-series data
% X0    Snapshot of Xf before last 'abs', 'pow', or 'env' operation
%
% SPECIFYING THE RULE SET
%   A cell string; rules and options are applied in sequential order, so
% e.g. {'demean'; 'env'} and {'env'; 'demean'} are NOT identical.
%
% RULE                      PROCEDURE
% 'demean'                  Remove mean of each column
% 'detrend'                 Detrend each column
% 'env'                     Compute envelope of each column
% 'filt',[corners],rule     Apply 2-pole butterworth zero-phase filter
%                               'filt' *must* be followed by:
%                               --corner frequencies in Hz (e.g. [1 12])
%                               --filtering rule (e.g. 'bandpass')
% 'mir'                   vertcat(X,flipud(X));
% 'abs'                   Absolute value
% 'pow',p                 Raise X to power p
%
% ========================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.0, 2016-02-29

Xf = X;
X0 = X;
for j = 1:1:numel(r)
    if ischar(r{j})
        switch lower(r{j})
            case 'demean'
                Xf = detrend(Xf,'constant');
            case 'detrend'
                Xf = detrend(Xf);
            case 'env'
                X0 = Xf;
                Xf = abs(hilbert(Xf));
            case 'filt'
                fwin = r{j+1};
                frule = r{j+2};
                if strcmpi(frule,'')
                    [b1,a1] = butter(2,fwin*2/fs);
                else
                    [b1,a1] = butter(2,fwin*2/fs,frule);
                end
                Xf = filtfilt(b1,a1,Xf);
            case 'mir'
                Xf = vertcat(Xf, flipud(Xf));
            case 'abs'
                X0 = Xf;
                Xf = abs(Xf);
            case 'pow'
                X0 = Xf;
                Xf = abs(Xf).^r{j+1};
            otherwise
        end
    end
end
