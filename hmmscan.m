function Ev = hmmscan(X, varargin)
% Ev = hmmscan(X, OPTS)
%
%       Scan data X for events using a 2-state HMM with an expectation-
% maximization algorithm to classify state memberships.
%
% INPUTS
% X     Data arranged in column vectors
% OPTS  Program options in 'name',value pairs, e.g. 'Lw', 3000. See below
%       for defaults.
% 
% OUTPUTS
% Ev    Event structure:
%           s   Event start time, in seconds (relative to scan begin)
%           e   Event end time
%           p   P-pick
%           pq  P-quality (see paper)
%           dp  P-pick error
%           rs  Smoothed detection statistic
%           rb  raw detection statistic
%
% OPTS      Default         Meaning
% Lw        2000            Long window, samples
% Sw        200             Short window, samples
% Th        1.3             Detection threshold (Sw >= Th*Lw) [1]
% cf        'env'           Characteristic function to transform X
% dc        1               Check for duplicates? (1 = 'yes')
% eos       500             Event offset, samples
% f0        0.10            Initial outlier fraction
% fs        4000            Sampling rate, Hz
% fw        [0.005 0.4]     Allowed range for outlier fraction
% minsep    25              Minimum offset of non-duplicate events, samples
% op        1               Pick onset? (1 = 'yes')
% pad       0               Pad start with esp, end with eep? 1 = 'yes'
%   eep     2000            Event end padding, samples
%   esp     2000            Event start padding, samples
% plot      0               Create plot? (1 = 'yes')
% slv       'exp'           HMM solver for expectation maximization
% s0        6               How many channels must converge to a solution?
% v         0               Verbosity
% weed      1               Weed traces whose HMM states are too similar?
%
% Dependencies: exmax.m, genprep.m, emtheta.m, es94ap.m, sop.m
% 
% NOTES
% [1]   See paper! Generally this should be set MUCH LOWER than for
% traditional STA-LTA; unless an environment is very noisy, Th > 1.5 is
% practically never needed.
% 
% ========================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.0, 2016-02-29

% Defaults for opts
opts.Lw = 2000;
opts.Sw = 200;
opts.Th = 1.3;
opts.cf = 'env';
opts.dc = 1;
opts.eep = 2000;
opts.eos = 1000;
opts.esp = 2000;
opts.f0 = 0.10;
opts.fs = 4000;
opts.fw = [0.005 0.4];
opts.minsep = 25;
opts.op = 1;
opts.pad = 0;
opts.plot = 0;
opts.slv = 'exp';
opts.s0 = 6;
opts.v = 0;
opts.weed = 1;


% ______________________________________________________________________
% Parse varargin
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

% ________________________________________________________________
% Computations begin
X = X./max(abs(X(:)));
Y = genprep(X, {opts.cf});

% Some needed variables next
[Ny,Nc] = size(Y);
L = opts.Lw;
sw2 = floor(opts.Sw/2);
L1 = ceil(L/2);

% Initialize counters
yc = 0;                      % Y counter
ec = 0;                      % Event counter

% Set initial values
if numel(opts.f0) == 1
    opts.f0 = repmat(opts.f0,[1 Nc]);
end 
opts.Q0 = hmminit(Y(1:L,:),opts.f0);

% Initialize structures
Ev.rs = zeros(size(Y));
Ev.rb = zeros(Ny,1);
Ev.thr = zeros(1,Nc);
Ev.pt = Ev.rs;

% Sort options structure
opts = orderfields(opts);

% Begin
if opts.v > 0
    tic;
    fprintf(1,'Working...%5.1f%s',0,'%');
end
while yc+L < Ny
    L = min(L,Ny-yc);
    y = Y(yc+1:yc+L,:);
    Ev.thr(yc+1:Ny,:) = repmat(opts.Th.*opts.f0,[Ny-yc 1]);
    
    % Evaluate HMM
    [pt,f,Q] = exmax(y,opts);
    nsol = numel(f);
    
    % Proceed if nsol > solmin
    if nsol > opts.s0
        key = Q{2};
        Ev.pt(yc+1:yc+L,key) = pt;
        Ev.thr(yc+1:Ny,key) = repmat(opts.Th.*f,[Ny-yc 1]);
        R = sop(pt,f,opts);
        
        % Check for events
        if ~isempty(R.s)
            J1 = numel(R.s);
            
            % Add start time, end time, quality to database
            for j1 = 1:J1
                Ev.s(ec+j1,:) = R.s(j1) + yc;
                Ev.e(ec+j1,:) = R.e(j1) + yc;
                Ev.eq(ec+j1,key) = R.eq(j1,:);
                
                % Check for picks if picking
                if opts.op
                    if numel(find(R.p))
                        
                        % Add pick, quality, error bars to data base
                        Ev.p(ec+j1,key) = R.p(j1,:) + yc*(R.p(j1,:)>0);
                        Ev.pq(ec+j1,key) = R.pq(j1,:);
                        Ev.dp(ec+j1,key) = R.dp(j1,:);
                    else
                        Ev.p(ec+j1,key) = 0;
                        Ev.pq(ec+j1,key) = 0;
                        Ev.dp(ec+j1,key) = 0;
                    end
                end
                
                % Add smoothed R
                if isfield(R,'rs')
                    Ev.rs(yc+1:yc+L,key) = R.rs;
                    Ev.rs(yc+L+1:Ny,key) = repmat(min(R.rs),[Ny-yc-L 1]);
                    Ev.rb(yc+1:yc+L,1) = R.rb;
                    Ev.rb(yc+L+1:Ny,1) = repmat(min(R.rb),[Ny-yc-L 1]);
                end
            end
            
            ec = ec + J1;
            yc = yc + max([L1,max(R.e)]);
        else
            yc = yc + L1;
        end
        
        clear R;
        opts.f0(key) = min([max([f opts.fw(1)]) opts.fw(2)]);
        opts.S0(:,key) = Q{1};
    else
        yc = yc + L1;
    end
    if opts.v > 0
        fprintf(1,'\b\b\b\b\b\b%5.1f%s',100*yc/Ny,'%');
    end
end
if opts.v > 0
    fprintf(1,'\b\b\b\b\b\bdone! Elapsed time %4.2f s\n',toc);
end

%% Check for duplicates
d1 = 0;
if isfield(Ev,'s')
    if opts.dc
        j1 = 1;
        J1 = numel(Ev.s);
        while j1 < J1;
            if abs(Ev.s(j1+1)-Ev.e(j1)) < opts.minsep
                Ev.e(j1) = Ev.e(j1+1);
                Ev.eq(j1,:) = max([Ev.eq(j1,:); Ev.eq(j1+1,:)]);
                
                % Remove next event
                Ev.s(j1+1,:) = [];
                Ev.e(j1+1,:) = [];
                
                if max(Ev.eq(j1,:) > 0)
                    Ev.eq(j1+1,:) = [];
                    
                    if opts.op
                        Ev.p(j1+1,:) = [];
                        Ev.pq(j1+1,:) = [];
                        Ev.dp(j1+1,:) = [];
                    end
                else
                    Ev.eq(j1,:) = Ev.eq(j1+1,:);
                    Ev.eq(j1+1,:) = [];
                    
                    if opts.op
                        Ev.p(j1,:) =  Ev.p(j1+1,:) ;
                        Ev.pq(j1,:) = Ev.pq(j1+1,:);
                        Ev.dp(j1,:) = Ev.dp(j1+1,:);
                        Ev.p(j1+1,:) = [];
                        Ev.pq(j1+1,:) = [];
                        Ev.dp(j1+1,:) = [];
                    end
                end
                
                J1 = numel(Ev.s);
                d1 = d1+1;
            else
                j1 = j1+1;
            end
        end
        if opts.v && d1 > 0
            warning(['Merged ' num2str(d1) ' duplicates']);
        end
    end
    if opts.v
        disp(['Found ' num2str(numel(Ev.s)) ' good events']);
    end
end

%% Shift rs to align with midpoint of each detection window
Ev.rs = [zeros(sw2+1,Nc); Ev.rs(1:Ny-sw2-1,:)];
Ev.rb = [zeros(sw2+1,1); Ev.rb(1:Ny-sw2-1,1)];
Ev.time = toc;

% Set event start, end, P times to be in seconds relative to record begin
if isfield(Ev,'s')
    Ev.s = Ev.s/opts.fs;
    Ev.e = Ev.e/opts.fs;
    
    if isfield(Ev,'p')
        Ev.p = Ev.p/opts.fs;
        Ev.dp = Ev.dp/opts.fs;
        if size(Ev.p,2) < Nc
            Ev.p(:,end+1:Nc) = 0;
        end
    end
    
    %% Plot detections
    if opts.plot
        t = 1/opts.fs:1/opts.fs:size(Y,1)/opts.fs;
        figure;
        hold on;
        for k = 1:1:Nc
            x = k + X(:,k)./(2*max(abs(X(:,k))));
            plot(t,x,'k-','linewidth',1);
        end
        y1 = 0.5;
        y2 = size(X,2)+0.5;
        ylim([y1 y2]);
        hold on;
        for k = 1:1:numel(Ev.s)
            plot(Ev.s(k)*[1 1],[y1 y2],'g-');
            plot(Ev.e(k)*[1 1],[y1 y2],'r-');
           if isfield(Ev,'p')
               for n = 1:1:size(Ev.p,1)
                    for c = 1:1:Nc
                        if Ev.p(n,c) > 0
                            plot(Ev.p(n,c)*[1 1],[c-0.5 c+0.5],'b-');
                        end
                    end
                end
           end
        end
        axis tight;
    end
end
