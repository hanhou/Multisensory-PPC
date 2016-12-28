function [alpha,fval,dfval,nfe,nge] = wolfe_search3(f,df,x,p, ...
    alpha0,c1,c2,f0,df0,trace)
% [ALPHA,FVAL,DFVAL,NFE,NGE] = ...
%	wolfe_search2(F,DF,X,P,ALPHA0,C1,C2,F0,DF0,TRACE)
% Performs search for point satisfying the two strong Wolfe conditions:
%
%	f(x+\alpha p) <= f(x) + c1\alpha p'.grad f(0)		(WC1)
%	|p'.grad f(x+\alpha p)| <= c2|p'.grad f'(x)|		(WC2b)
%
% ALPHA0 is the initial value of alpha in the search.  If ALPHA0
%	satisfies the Wolfe conditions, then the search stops there. 
% F0 is f(x) at the initial x.
% DF0 is grad f(x) at the initial x.
% DF(PARAMS,X) is the gradient of F(PARAMS,X) with
%	respect to X
% Note that PARAMS is the parameter for F: F(PARAMS,X)
% If TRACE is non-zero then print out information about the process
%
% The algorithm follows that of Wright & Nocedal "Numerical Optimization",
% pp. 60-62, section 3.4
%
% Returns
%	ALPHA	-- step length parameter
%	FVAL	-- function value at end
%	DFVAL	-- gradient value at end
%	NFE	-- # function evaluations
%	NGE	-- # gradient evaluations
%	modif by F KLAM/ functions without additional params

%
% Check parameters
%
alpha = 0;
fval = f0;
dfval = df0;
nfe = 0;
nge = 0;

if ~ ( 0 < c1 & c1 < c2 & c2 < 1 )
  fprintf('wolfe_search2: Error: Need 0 < c1 < c2 < 1\n');
  if trace ~= 0
    fprintf('wolfe_search2: c1 = %g, c2 = %g\n', c1, c2);
  end
  return
end

% Initialization
if nargin <= 10
  trace = 0;
end

slope0 = p'*df0;
if slope0 >= 0
  fprintf('wolfe_search2: Error: Need a descent direction\n');
  if trace ~= 0
    fprintf('wolfe_search2: p''.grad f(x) = %g\n',slope0);
  end
  return
end

%
% Bracketing phase
%
if trace ~= 0
  fprintf('wolfe_search2: Bracketing phase: alpha = %g\n', alpha0);
end
alpha = alpha0;
old_alpha = 0;
% f0  = feval(f, params,x);
% df0 = feval(df,params,x);
old_fval = f0;
old_dfval = df0;
old_slope = slope0;
if trace ~= 0
  fprintf('wolfe_search2: f(x) = %g, p''.grad f(x) = %g\n',f0,slope0);
end

% Main loop
firsttrip = 1;
while 1 %%% forever do...
  xplus = x+alpha*p;
  fval = feval(f,xplus); %fval = feval(f,params,xplus);
  nfe = nfe + 1;
  if trace ~= 0
    fprintf('wolfe_search2: alpha = %g, f(x+alpha*p) = %g\n', alpha, ...
	fval);
  end
  if ( fval > f0 + c1*alpha*slope0 ) | ( ( ~ firsttrip ) & ...
	( fval >= old_fval ) ) 
    if trace ~= 0
      fprintf('wolfe_search2: (WC1) failed or f increased\n');
    end
    break;
  end
  if trace ~= 0
    fprintf('wolfe_search2: (WC1) holds & f decreased\n');
  end
  dfval = feval(df,xplus); %dfval = feval(df,params,xplus);
  nge = nge + 1;
  slope = p'*dfval;
  if trace ~= 0
    fprintf('wolfe_search2: p''.grad f(x+alpha*p) = %g\n', slope);
  end
  if ( abs(slope) <= c2*abs(slope0) )
    if trace ~= 0
      fprintf('wolfe_search2: (WC2) holds\n');
    end
    return
  end
  if ( slope >= 0 )
    if trace ~= 0
      fprintf('wolfe_search2: f''(alpha) >= 0\n');
    end
    break;
  end

  % Update variables -- note no upper limit on alpha
  temp = alpha;
  % alpha = alpha + old_alpha;
  alpha = 2*alpha;
  old_alpha = temp;
  old_fval = fval;
  old_slope = slope;
end

if ( trace ~= 0 )
  fprintf('wolfe_search2: Entering ''zoom'' phase.\n')
end

%
% "Zoom" phase
%
alpha_lo = old_alpha;
alpha_hi =     alpha;
f_lo     = old_fval;
f_hi     =     fval;
df_lo    = old_dfval;
slope_lo = old_slope;
%df_hi   =     dfval;
%slope_hi=     slope;

if trace ~= 0
  fprintf('wolfe_search2: zoom phase: alpha_lo = %g, alpha_hi = %g\n', ...
      alpha_lo, alpha_hi);
end

iter_cnt = 0;
while 1 %%% forever do...

  % form quadratic interpolant of function values on alpha_lo, alpha_hi
  % and the derivative at alpha_lo...
  % and find min of interpolant within interval [alpha_lo,alpha_hi]
  a = f_lo;
  b = slope_lo;
  dalpha = alpha_hi-alpha_lo;
%if dalpha^2==0
%  	return;
%end
  c = (f_hi - f_lo - dalpha*slope_lo)/dalpha^2;
  % iter_cnt
  if ( ( c <= 0 ) | ( mod(iter_cnt,3) == 2 ) )
    % Use bisection
    alpha = alpha_lo + 0.5*dalpha;
    if trace ~= 0
      fprintf('wolfe_search2: using bisection: c=%g\n', c);
    end
  else
    % Use min of quadratic
    alpha = alpha_lo - 0.5*b/c;
    if trace ~= 0
      fprintf('wolfe_search2: using quadratic: a=%g, b=%g, c=%g\n', ...
	  a, b, c);
    end
  end
  if trace ~= 0
    fprintf('wolfe_search2: alpha = %g\n', alpha);
  end

  % main part of loop
  xplus = x + alpha*p;
  fval = feval(f,xplus); %fval = feval(f,params,xplus);
  nfe = nfe + 1;
  if trace ~= 0
    fprintf('wolfe_search2: f(x+alpha*p) = %g\n', fval);
  end
  if ( ( fval > f0 + c1*alpha*slope0 ) | ( fval >= f_lo ) )
    if trace ~= 0
      fprintf('wolfe_search2: zoom: (WC1) fails or f increased\n');
      fprintf('wolfe_search2: zoom: update alpha_hi\n');
    end
    alpha_hi = alpha;
    f_hi   = fval;
    % df_hi = feval(df,params,xplus);
  else
    if trace ~= 0
      fprintf('wolfe_search2: zoom: (WC1) holds and f decreased\n');
    end
    dfval = feval(df,xplus); %dfval = feval(df,params,xplus);
    nge = nge + 1;
    slope = p'*dfval;
    if trace ~= 0
      fprintf('wolfe_search2: p''.grad f(x+alpha*p) = %g\n', slope);
    end
    if ( abs(slope) <= c2*abs(slope0) )
      if trace ~= 0
	fprintf('wolfe_search2: zoom: (WC2b) holds\n');
      end
      return
    end
    if ( slope*dalpha >= 0 )
      if trace ~= 0
	fprintf('wolfe_search2: zoom: alpha_hi <- alpha_lo & update alpha_lo\n');
      end
      alpha_hi = alpha_lo;
      f_hi     = f_lo;
      % df_hi  = df_lo;
      % slope_hi= slope_lo;
      alpha_lo = alpha;
      f_lo     = fval;
      df_lo    = dfval;
      slope_lo = slope;
    else
      if trace ~= 0
	fprintf('wolfe_search2: zoom: update alpha_lo\n');
      end
      alpha_lo = alpha;
      f_lo     = fval;
      df_lo    = dfval;
      slope_lo = slope;
    end
  end
  % Update iteration count
  iter_cnt = iter_cnt + 1;

end % of forever do...

