function [x,fval,g,nfe,nge,outflag] = ...
    cg_pr(func,dfunc,x,tol,ftol,c1,c2,it_limit,trace)
% [X,FVAL,GVAL,NFE,NGE] = cg_pr(FUNC,DFUNC,PARAMS,X,TOL,C1,C2,IT_LIMIT,TRACE)
% Polak-Ribiere conjugate gradient method (without restart)
% Uses Wolfe line-search conditions (parameters c1 and c2)
%
% FUNC is the name of the function to minimize: FUNC(PARAMS,X)
% DFUNC is the name of the gradient function: DFUNC(PARAMS,X)
% PARAMS are parameters that are passed to FUNC and DFUNC
% X is the current position
% TOL is the stopping tolerance: ||gradient|| < tol
% C1 and C2 are the Wolfe line-search parameters
% IT_LIMIT is the iteration limit
% If TRACE is non-zero, then print out information about process
%
% Returns
%	X	-- final "minimizer"
%	FVAL	-- FUNC(PARAMS,X), function value at X
%	GVAL	-- DFUNC(PARAMS,X), gradient value at X
%	NFE	-- # function eval'ns
%	NGE	-- # gradient eval'ns
%%%modified by FKLAM outflag / functions without arguments other than x...

% Initialization
if nargin < 8
  trace = 0;
end

n = length(x);
pfval=1000;

fval = feval(func,x); %fval = feval(func,params,x);
g = feval(dfunc,x); %g = feval(dfunc,params,x);
nfe = 1;
nge = 1;
p = -g;
g_sqr = g'*g;
for k = 0:it_limit
  [alpha,fval,new_g,nfe_ws,nge_ws] = ...
      wolfe_search3(func,dfunc,x,p,1.0,c1,c2,fval,g,trace);
  nfe = nfe + nfe_ws;
  nge = nge + nge_ws;
  if trace ~= 0
    fprintf('cg_fr: iter # %d, alpha = %g, function value = %g\n', ...
	k, alpha, fval);
  end
  x = x + alpha*p;
  new_g_sqr = new_g'*new_g;
  % Polak Ribiere modification
  beta = new_g'*(new_g-g)/g_sqr;
  g_sqr = new_g_sqr;
  g = new_g;
  if norm(g,2) < tol | abs(pfval-fval)<ftol
  outflag=(k<it_limit);
    return
  end
  % Restart code (if desired)
  % if k % n == n-1
  %   beta = 0;
  % end
  p = -g + beta*p;
  pfval=fval;
end
%fprintf('Have used %d iterations; returning without finding solution\n', ...
%    it_limit);
outflag=(k<it_limit);

