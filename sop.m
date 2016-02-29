function R = sop(p, f, varargin)
% R = sop(P, F, OPTS);
%
% Compute the short-term-averaged outlier probability given outlier
% probability time series P and outlier fraction F. If P is a matrix and F
% is an array, it's required that numel(F) = size(P,2). 
% 
% INPUTS
% P         outlier probability time series, arranged as column vectors
% F         outlier fraction for each column of P
% 
% OUTPUT STRUCTURE R
% FIELD         MEANING
%   rb          Detection statistic
%   rs          Smoothed rho_t from onset picking
%   t           Central sample # of each rho_t calculation
%   s           Event start times, samples
%   e           Event end times, samples
%   p           Event P picks, samples

%
% OPTS  DEFAULT         MEANING
% Sw    100             Short window, samples
% Th    1.5             Event windows correspond to STA - Th*LTA > 0
% pTh   2.0             Only autopick if STA - pTh*LTA > 0
% eos   2000            Event start times must be > eos apart
% op    0               Do onset picking (1 = 'yes')
% pad   0               Pad event start with esp, end with eep? 1 = 'yes'
%  esp  2000            Event start pad, in samples
%  eep  4000            Event end pad, in samples
% qth   0.1             Only retain picks with q > qth
% tmax  0.95            Events with thr > tmax are set to tmax
% v     0               Verbosity
%
% ========================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.0, 2016-02-29


%% Default options
opts.Sw = 200;
opts.Th = 1.3;
opts.pTh = 2.0;
opts.eep = 2000;         % Event end (post) pad, in samples
opts.eos = 1000;         % Event start times must be > eos apart
opts.esp = 2000;         % Event start (pre) pad, in samples
opts.op = 0;             % Onset picking, 1 = 'yes'
opts.pad = 0;            % Pad start with esp, end with eep? 1 = 'yes'
opts.qth = 0.1;          % Only try to pick events with q > qth
opts.tmax = 0.95;        % thr > opts.tmax are set to opts.tmax
opts.v = 0;              % Verbosity

% Dependent defaults
opts.dt = opts.Sw/2;     % Time between STA calculations
opts.dpmax = opts.Sw;    % P picks with errors over dpmax aren't saved

%% Parse varargin
if numel(varargin) > 0
    j = 1;
    while j <= numel(varargin)
        if isstruct(varargin{j})
            inputs = fieldnames(varargin{j});
            for k = 1:numel(inputs)
                eval(['opts.' inputs{k} '= varargin{j}.' ...
                      inputs{k} ';']);
            end
            j = j+1;
        else
            eval(['opts.' varargin{j} '= varargin{j+1};']);
            j = j+2;
        end        
    end
end
opts = orderfields(opts);

%% Main
% Create vector t of times at which we compute sop
[Np,Nc] = size(p);
t = 0:opts.dt:Np-opts.Sw;
K = numel(t);
sw2 = ceil(opts.Sw/2);

% Set event detection options if we are doing event detection
s1 = [];
e1 = [];
eq = [];
evip = 0;
m = 0;
eos = opts.eos;

% Ensure thresholds are correctly specified for each channel
if numel(opts.Th) < Nc
    thr = opts.Th;
    thr(end+1:Nc) = repmat(thr(end),[1 Nc-numel(opts.Th)]);
    if opts.op
        pthr = opts.pTh;
        pthr(end+1:Nc) = repmat(thr(end),[1 Nc-numel(opts.pTh)]);
    end
else
    thr = opts.Th;
    if opts.op
        pthr = opts.pTh;
    end
end

% thr<1 is an exact probability p0; thr >= 1 is a scalar multiplier of f
if min(thr) >= 1
    thr = thr.*f;
    if opts.op
        pthr = pthr.*f;
    end
end

% Ensure thr are sane
thr(thr>opts.tmax) = opts.tmax;           
if opts.op
    pthr(pthr>opts.tmax) = opts.tmax;
end

%% Compute detection statistic r
for k = 1:1:K
    rt = p(t(k)+1:t(k)+opts.Sw,:);
    r(k,:) = mean(rt)-thr;
end

%% Event detection
rs = zeros(size(p));
rb = zeros(Np,1);

% Initiatlize picking variables if doing onset picks
if opts.op
    op = zeros(1,Nc);
    pq = zeros(1,Nc);
    dp = zeros(1,Nc);
end

% Event detection at each time in array t
for k = 1:1:K
    % Fill rs
    rs(t(k)+1:t(k)+opts.Sw,:) = repmat(r(k,:),[opts.Sw 1]);
    
    % Event tracker
    etr = mean(r(k,:));
    
    % Stored for plotting purposes
    rb(t(k)+1:t(k)+opts.Sw,1) = etr;
    
    if etr > 0 & ~evip
        evip = 1;
        if ~m | t(k)-s1(m) >= eos
            m = m+1;
            s1(m,1) = t(k);
            k0 = k;
        end
    elseif (etr <= 0 | k == K) & evip
        evip = 0;
        
        % Compute event quality
        eq(m,:) = max(max(r(k0:k,:))./(1-thr),0);
        
        % Onset picking if selected & feasible. The condition numel(e1) <
        % m prevents repicking when end times are extended.
        if opts.op & numel(e1) < m
            op(m,:) = zeros(1,Nc);
            pq(m,:) = zeros(1,Nc);
            dp(m,:) = zeros(1,Nc);
            
            % ensure we aren't picking in coda of previous event
            if m > 1
                tt = e1(m-1)+opts.eep;
            else
                tt = 0;
            end
            if s1(m,1)-sw2 > tt & t(k)+sw2 < Np
                ps = s1(m,1)-sw2;
                pe = t(k)+sw2;
                p0 = p(ps:pe,:);
                [op(m,:),pq(m,:),dp(m,:),rs(ps+sw2:pe-sw2,:)] = ...
                    autopick(p0,pthr,ps,opts);
            end
        end
        e1(m,1) = t(k);
        
    end
end

% Check for event start at last time
if numel(s1) > numel(e1)
    eq(m,:) = max(r(k0:k,:))./(1-thr);
    if opts.op & numel(e1) < m
        if m > 1
            %tt = e1(m-1)+opts.eos*opts.eep;
            tt = e1(m-1)+opts.eep;
        else
            tt = 0;
        end
        if s1(m,1)-sw2 > tt & t(k)+sw2 < Np
            ps = s1(m,1)-sw2;
            pe = t(k)+sw2;
            [op(m,:),pq(m,:),dp(m,:),rs(ps+sw2:pe-sw2,:)] = ...
                autopick(p(ps:pe,:),thr,ps,opts);
        end
    end
    e1(m,1) = s1(m,1)+opts.Sw;
end

% Adjust event start, times
if opts.pad
    s1 = s1 - opts.esp;   % Event start pad
    e1 = e1 + opts.eep;   % Event end pad
end

s1 = s1 + round(opts.Sw/2) + 1;
e1 = e1 + round(opts.Sw/2) + 1;

% Fill rest of rs, rb
rs(t(K)+opts.Sw+1:Np,:) = repmat(min(r),[Np-t(K)-opts.Sw,1]);
rb(t(K)+opts.Sw+1:Np,1) = repmat(min(rb),[Np-t(K)-opts.Sw,1]);

%% Fill results structure
R = opts;
R.t = t + sw2;
R.r = r; 
R.s = s1;
R.e = e1;
R.eq = eq;
R.evip = evip;
R.rs = rs;
R.rb = rb;

% Control for poor quality picks
if opts.op & numel(s1)
    Nev = numel(R.s);
    op(dp > opts.dpmax) = 0;
    pq(dp > opts.dpmax) = 0;
    dp(dp > opts.dpmax) = 0;
    [Np,Ncp] = size(op);
    R.p = op;
    R.pq = pq;
    R.dp = dp;
    if Np < Nev | Ncp < Nc
        R.p(Np+1:Nev,Ncp+1:Nc) = 0;
        R.pq(Np+1:Nev,Ncp+1:Nc) = 0;
        R.dp(Np+1:Nev,Ncp+1:Nc) = 0;
    end
end

% ========================================================================
function [op,pq,dp,xs] = autopick(p0,thr,ps,opts)
[np,nc] = size(p0);
op = zeros(1,nc);
pq = zeros(1,nc);
dp = zeros(1,nc);
sw2 = round(opts.Sw/2);

for k1 = 1:1:np-opts.Sw
    pt(k1,:) = mean(p0(k1+1:k1+opts.Sw,:)) - thr;
end
[tp,oq,xs] = es94ap(pt);  % Autopick
oq = oq./(1-thr);         % Max possible value of 1

for c = 1:nc
    if tp(c) > 0 & oq(c) >= opts.qth
        op(c) = 1 + ps + tp(c) + sw2;       % Time correction
        pq(c) = oq(c);                      % Pick quality
        dp(c) = max(7,exp(1/pq(c)));        % Pick error
    else
        op(c) = 0;
        pq(c) = 0;
        dp(c) = 0;
    end
end
