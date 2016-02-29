function Q = hmminit(X, f, varargin)
% Q = hmminit(X, f, sord);
%
%   Generate initial values for the controlling parameters Q (mean and
% variance) of a 2-state HMM; the outlier population uses the outermost f%
% of X values.
%
% INPUTS
% X     Preprocessed data in column vectors
% f0    Initial outlier fraction
% sord  Sort order
% 
% OUTPUTS
% t     Parameters theta for test distributions, with row ordering:
%       1   m1      Mean of outlier population
%       2   m0      Mean of null population
%       3   v1      Variance of outlier population
%       4   v0      Variance of null population
%
% ========================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.0, 2016-02-29

s = {'abs'; 'descend'};
if numel(varargin) > 0
    if ~isempty(varargin{1})
        s = varargin{1};
    end
end
if numel(s) == 1
    s{2} = 'descend';
end

Nx = size(X,1);
F0 = round(f*Nx);

% Generate trial values using a sorted X
switch lower(s{1})
  case 'abs'
    X0 = sort(abs(X),1,s{2});
  case 'sq'
    X0 = sort(X.^2,1,s{2});
  otherwise
    X0 = sort(X,1,s{2});
end
    
Q = [mean(X0(1:F0,:)); mean(X0(F0+1:end,:)); ...
     var(X0(1:F0,:)); var(X0(F0+1:end,:))];