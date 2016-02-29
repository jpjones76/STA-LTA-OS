function varargout = exmax(X, varargin)
% P = exmax(X);
% P = exmax(X,OPTIONS);
% [P,f,Q,K] = exmax(X,OPTIONS);
%
% Two-state expectation-maximization for classification. Compute the
% probability P that each point in X belongs to one of two HMM states
% for each channel, having outlier fraction f and state parameters Q.
%
% If X is a matrix, expectation-maximization is applied to each column
% independently. 
%
% OUTPUTS -- P is required, rest optional
% P     Probability that each point in X is an outlier.
% f     Fraction belonging to outlier population in each column.
% Q     Controlling statistics (mean and variance, if applicable)
%       for each population of each column. 
%       S[1:2] = mean (1=null, 2=outlier). S[3:4] = var. (same order)
% K     Key giving correspondence between columns in P, f, Q and columns
%       in the input
%
% OPTIONS
%   Specify options with 'name', value pairs, e.g. 'f0',0.1.
%
% Option    Default     Meaning
% Ni        100         Number of iterations
% f0        0.1         Initial outlier fraction
% slv       'exp'       Solver to use [1]
% Q0        [unset]     Initial values for controlling parameters Q
% tol       -1          Tolerance for convergence (-1 = set adaptively)
% v         0           Verbosity
% weed      1           Weed solutions with identical states? 1 = "yes"
%   nstd    3           If weeding, states within nstd are identical
%   cf      'abs'       Characteristic function that transformed data [1]
% 
% [1] cf and slv are cross-referenced in MinThr.mat; a bad cf string will
% only affect solution weeding.
% 
% ========================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.0, 2016-02-29

[Nx,Nc] = size(X);

% Default options
opts.Ni = 100;
opts.cf = 'env';
opts.f0 = 0.10;
opts.nstd = 3;
opts.slv = 'exp';
opts.th = 0; 
opts.tol = -1;
opts.v = 0;
opts.weed = 1;

%% Initialization
% Parse varargin
if nargin > 0
    j = 1;
    while j <= numel(varargin)
        if isstruct(varargin{j})
            inputs = fieldnames(varargin{j});
            for k = 1:numel(inputs)
                eval(['opts.' inputs{k} '= varargin{j}.' inputs{k} ';']);
            end
            j = j+1;
        else
            eval(['opts.' varargin{j} '= varargin{j+1};']);
            j = j+2;
        end        
    end
end
opts = orderfields(opts);

% Set verbosity, start timing if v > 0
v = opts.v;
if v > 1; tic; end

% Set tolerance dynamically based on single precision
if opts.tol == -1
    tol = 1.0e-6*Nx;
else
    tol = opts.tol;
end

if numel(opts.f0) < Nc
    f = repmat(opts.f0(1),[1 Nc]);
else
    f = opts.f0;
end

% Set initial outlier fractions, if needed
if isfield(opts,'Q0')
    Q0 = opts.Q0;
else
    if ~isfield(opts,'sort')
        if v>1; warning('opts.sort not found, using default'); end
        switch lower(opts.slv)
            case 'beta'
                opts.sort = {'raw'; 'ascend'};
            otherwise
                opts.sort = {'abs'; 'descend'};
        end
    end
    Q0 = hmminit(X, opts.f0, opts.sort);
    if strcmpi(opts.slv,'ray')
        Q0(1:2,:) = (Q0(1:2,:).^2)./pi;
    end
end
m1 = Q0(1,:);           % Mean of pop 1
m0 = Q0(2,:);           % Mean of pop 0
v1 = Q0(3,:);           % Variance of pop 1
v0 = Q0(4,:);           % Variance of pop 0
i1 = 1:1:Nc;            % Indices to columns
mm = 0;                 % Index to last solution

Af = [];
ff = [];
m1f = []; 
m0f = [];
v1f = []; 
v0f = [];
key = [];

%% ExMax iterative classifier
L0 = zeros(1,Nc);
for nn = 1:1:opts.Ni
    [L,A] = exstep(opts.slv,X,f,m1,m0,v1,v0);
    dL = L-L0;
        
    % Pull solutions before we proceed
    if nn > 1
        ind = find(dL<tol); 
        if ~isempty(ind)
            if v > 1
                disp(['At n=' num2str(nn) ', these chans. converged:']);
                disp(num2str(i1(ind)));
            end
            for k = 1:numel(ind)
                k1 = ind(k);
                mm = mm+1;
                key(mm) = i1(k1);
                
                % Store information in variables denoted "f" for "final"

                Af(:,mm) = A(:,k1); ff(mm) = f(k1);
                m1f(mm) = m1(k1); m0f(mm) = m0(k1);
                v1f(mm) = v1(k1); v0f(mm) = v0(k1);
            end
            
            % Remove solved channels
            A(:,ind)= []; 
            X(:,ind)= [];
            L(ind)  = []; 
            i1(ind) = [];
        end
    
        % If all solutions are found, or we're at nn = Ni, we are done.
        if (isempty(L) || sum(isnan(L)) == numel(L)) && v > 1
            disp('All solutions found! Order:');
            disp(num2str(key));
            break
        elseif nn == opts.Ni && v > 1 && ~isempty(L)
            Warning(['Some channels failed to converge after '...
                  num2str(opts.Ni) ' iterations. Exiting.']);
            if v > 1
                disp('Failures:');
                disp(num2str(i1));
            end
            break
        end
    end
   
    % Maximization step: Set values for next iteration    
    f = sum(A)/Nx;
    if strcmpi(opts.slv,'ray')
        m1 = 0.5*sum(A.*X.^2)./sum(A);
        m0 = 0.5*sum((1-A).*X.^2)./sum(1-A);
    else
        m1 = sum(A.*X)./sum(A);
        m0 = sum((1-A).*X)./sum(1-A);
    end
    v1 = sum(A.*(X-repmat(m1,[Nx 1])).^2)./sum(A);
    v0 = sum((1-A).*(X-repmat(m0,[Nx 1])).^2)./sum(1-A);
    
    L0 = L;

    if v > 2
        disp(['Done iteration ' num2str(nn)]);
    end
end

%% Output control
% Output only channels that converged
if isempty(key)
    A = [];
    f = [];
    Q{1} = {};
    Q{2} = Q{1};
else
    [~,i2] = sort(key);
    A = Af(:,i2);
    f = ff(i2);
    Q{1} = [m1f(:,i2); m0f(:,i2); v1f(:,i2); v0f(:,i2)];
    Q{2} = key(i2);
end
clear Af ff Xf x0;

% Convert beta distribution mean, variance to parameters a, b
if strcmpi(opts.slv,'beta')
    m1 = Q{1}(1,:);
    m0 = Q{1}(2,:);
    v1 = Q{1}(3,:);
    v0 = Q{1}(4,:);    
    a1 = m1.*(((m1.*(1-m1))./v1)-1);
    b1 = (1-m1).*(((m1.*(1-m1))./v1)-1);
    a0 = m0.*(((m0.*(1-m0))./v0)-1);    
    b0 = (1-m0).*(((m0.*(1-m0))./v0)-1);
    Q{3} = [a0; b0; a1; b1];
end

%% Weed solutions whose null and outlier states are too similar
if opts.weed && ~isempty(f)
    load MinThr.mat
    [J1,K1] = size(MinThr.thmin);
    for j1 = 1:J1
        for k1 = 1:K1
            if (strcmpi(MinThr.pdf{j1,k1},opts.slv) * ...
                    strcmpi(MinThr.cf{j1,k1},opts.cf))
                t0 = MinThr.thmin(j1,k1) + opts.nstd*MinThr.thstd(j1,k1);
            end
        end
    end
    [th1,th2] = em_get_th(Q{1},f);
    switch lower(opts.slv)
        case {'beta','gauss'}
            th = 0.5 * (th1+th2); % Not recommended
        case 'zgauss'
            th = th2;
        otherwise
            th = th1;
    end
    ch = find(th > t0);

    Q{1} = Q{1}(:,ch);
    Q{2} = Q{2}(:,ch);
    A = A(:,ch);
    f = f(ch);
elseif isempty(f)
    warning('ExMax did not converge for any channel; weeding skipped!');
end

%% Parse varargout
varargout{1} = A;
if nargout > 1
    varargout{2} = f;
    if nargout > 2
        varargout{3} = Q;
    end
end

% End timing
if v > 1
    t1 = toc;
    disp(['exmax.m finished after ' num2str(t1) ' seconds.']);
end

% ========================================================================
% End main routine

% ========================================================================
%% Expectation step
function [L,A] = exstep(slv,X,f,m1,m0,v1,v0)
Nx = size(X,1);

switch lower(slv)
  case 'beta'
    b1 = (v1.*m1 - v1 + m1.^3 - 2*m1.^2 + m1)./v1;
    a1 = b1.*m1./(1-m1);
    b0 = (v0.*m0 - v0 + m0.^3 - 2*m0.^2 + m0)./v0;
    a0 = b0.*m0./(1-m0);    
    p1 = X.^(repmat(a1-1,[Nx 1])) ...
         .* (1-X).^(repmat(b1-1,[Nx 1])) ...
         ./ repmat(beta(a1,b1),[Nx 1]);
    p0 = X.^(repmat(a0-1,[Nx 1])) ...
         .* (1-X).^(repmat(b0-1,[Nx 1])) ...
         ./ repmat(beta(a0,b0),[Nx 1]);
    A = getA(p1,p0,f);
    L = sum(A).*log(f./beta(a1,b1)) ...
        + (a1-1).*sum(A.*log(X)) ...
        + (b1-1).*sum(A.*log(1-X)) ...
        + sum(1-A).*log((1-f)./beta(a0,b0)) ...
        + (a0-1).*sum((1-A).*log(X)) ...
        + (b0-1).*sum((1-A).*log(1-X));
    
  case 'ray'
    s1 = 1./m1;     % Inverse scale parameters speed up computations
    s0 = 1./m0;     
    p1 = X.*repmat(s1,[Nx 1]).* ...
         exp(-0.5.*repmat(s1,[Nx 1]).*X.^2);
    p0 = X.*repmat(s0,[Nx 1]).* ...
         exp(-0.5.*repmat(s0,[Nx 1]).*X.^2);
    A = getA(p1,p0,f);
    L = sum(A).*log(f.*s1) + sum(A.*log(X)) ...
        - 0.5*s1.*sum(A.*(X.^2)) ...
        + sum(1-A).*log((1-f).*s0) + sum((1-A).*log(X)) ...
        - 0.5*s0.*sum((1-A).*(X.^2));
    
  case 'zgauss'
    s1 = 1./sqrt(2.*v1);
    s0 = 1./sqrt(2.*v0);
    p1 = repmat(s1./sqrt(pi),[Nx 1]).* ...
         exp(-X.^2.*repmat(s1.^2,[Nx 1]));
    p0 = repmat(s0./sqrt(pi),[Nx 1]).* ...
         exp(-X.^2.*repmat(s0.^2,[Nx 1]));
    A = getA(p1,p0,f);
    L = sum(A).*(log(s1.*f./sqrt(pi))) + ...
        sum(1-A).*(log(s0.*(1-f)./sqrt(pi))) - ...
        s1.^2.*sum(A.*X.^2) - ...
        s0.^2.*sum((1-A).*X.^2);
    
  case 'gauss'
    s1 = 1./sqrt(2.*v1);
    s0 = 1./sqrt(2.*v0);
    p1 = repmat(s1./sqrt(pi),[Nx 1]).* ...
         exp(-(X-repmat(m1,[Nx 1])).^2.*repmat(s1.^2,[Nx 1]));
    p0 = repmat(s0./sqrt(pi),[Nx 1]).* ...
         exp(-(X-repmat(m0,[Nx 1])).^2.*repmat(s0.^2,[Nx 1]));
    A = getA(p1,p0,f);
    L = sum(A).*(log(s1.*f./sqrt(pi))) + ...
        sum(1-A).*(log(s0.*(1-f)./sqrt(pi))) - ...
        s1.^2.*sum(A.*(X-repmat(m1,[Nx 1])).^2) - ...
        s0.^2.*sum(A.*(X-repmat(m1,[Nx 1])).^2);
    
  otherwise
    % Default case: Exponential distributions
    im1 = 1./m1;
    im0 = 1./m0;
    p1 = exp(-X.*repmat(im1,[Nx 1])) .* ...
         repmat(im1,[Nx 1]);
    p0 = exp(-X.*repmat(im0,[Nx 1])) .* ...
         repmat(im0,[Nx 1]);    
    A = getA(p1,p0,f);
    L = log(f.*im1).*sum(A) ...
        - sum(A.*X).*im1 ...
        + log((1-f).*im0).*sum(1-A) ...
        - sum((1-A).*X).*im0;
end

% ========================================================================
function A = getA(p1,p0,f)
[Nx,~] = size(p1);
A = repmat(f,[Nx 1]) .* p1 ...
    ./ (repmat(f,[Nx 1]).*p1 + repmat(1-f,[Nx 1]).*p0);

% ========================================================================
function [th1,th2] = em_get_th(M,f)
m1 = f.*M(1,:) + (1-f).*M(2,:);
m2 = f.*M(3,:) + (1-f).*M(4,:);
th1 = sqrt( (((M(1,:)-m1).^2 + (M(2,:)-m1).^2)) ./ m1.^2);
th2 = sqrt( (((M(3,:)-m1).^2 + (M(4,:)-m1).^2)) ./ m2.^2);

