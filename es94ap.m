function [p,q,Rs] = es94ap(R, varargin)
% p = es94ap(R);
%
%    Autopick subtractive "ratio" time-series R with method modified from
% Earle and Shearer (1994). For multichannel data, R can be a matrix with
% per-channel values in column vectors.
%    The first pick of each channel k is returned as p(k). 
%
% p = es94ap(R,nh);
%    Specify hanning window length, nh, in samples. Default = 30.
% 
% [p,q] = es94ap(R);
%    Also return a quality q for the pick, equal to the value of the first
% local maximum. 
%
% [p,q,Rs] = es94ap(R);
%    Also return the smoothed ratio time series Rs.
% 
% Description
% 1.) At each channel k, smooth (lowpass) R(:,k) by filtering with a
% Hanning window. 
% 2.) If max(Rs) > 0, determine first local maximum tm1 in Rs. 
% 3.) From the first positive value of Rs to xm1, search for the last
% inflection point in Rs *before* tm1. This is the pick.
% 
% Non-Canonical Behavior
%    Unlike the original algorithm, we look first for upward inflection
% points, indicating sharp increases in Rs. 
%
% ========================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.0, 2016-02-29

% Defaults
nh = 30;               % Length of Hanning window

% Parse varargin
if numel(varargin) > 0
    if ~isempty(varargin{1})
        nh = varargin{1};
    end
end
K = size(R,2);
Rs = hannfilt(R, nh);

% Initialize p, q to always return a value, even a null one.
p = zeros(1,K);
q = p;

% Channel loop begins
for k = 1:1:K
    p(k) = 0;
    
    % Find t1 == first point where smoothed r exceeds 0. 
    t1 = find(Rs(:,k)>0,1,'first');
    if ~isempty(t1)
        
        % Take diff; time now begins at t1+1
        rd1 = diff(Rs(t1:end,k));
        
        % Find tm1 == first maximum in xs
        tm1 = find(rd1(2:end)<0 & rd1(1:end-1)>0,1,'first');
               
        % Does a maximum exist?
        if ~isempty(tm1)
            % Save as a measure of quality
            q(k) = Rs(tm1+1,k);
            
            % From t1 to tm1, compute xd2 = diff(xd1(t1:tm1));
            rd2 = diff(rd1(1:tm1));
            
            % Get last upward inflection point before tm1
            [~,p1] = find(rd2(2:end)<0 & rd2(1:end-1)>0,1,'last');
            if isempty(p1)
                [~,p1] = find(rd2(2:end)>0 & rd2(1:end-1)<0,1,'last');
            end            
            
            if numel(p1)
                p(k) = 2 + t1 + p1;
            end
        end
    end
    clear t1 tm1;
end


% =======================================================================
function xs = hannfilt(x, nh)
% Xs = hannfilt(X, nh);
%
%   Smooth vector or matrix of column vectors X by convolving with length
% nh hanning filter. nh is optional; default is 30 samples.

nx = size(x,1);
hf = 0.5*(1-cos(2*pi.*(0:1:nh-1)/(nh-1)));
hf = hf./sum(hf);
hfo = numel(hf)-1;
if nx <= 3*hfo
    x(nx+1:1+3*hfo,:) = 0;
end
xs = filtfilt(hf,1,x);
xs = xs(1:nx,:);
