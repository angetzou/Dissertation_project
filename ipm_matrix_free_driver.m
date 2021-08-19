function [x,y,s, optimal_solution] = ipm_matrix_free_driver(c,A,b,Q,k,maxit_ipm,maxit_kryl,printlevel,gamma,delta)
%IPDRIVER   for the matrix free Primal-dual interior-point method.
%
%  This is the driver function of a matrix free IPM for solving the
%   problems
%
%  min c'x + 0.5 x'Qx +0.5(x-xo)'Rp(x-xo) subject to A*x=b, x>=0,            (1)
%
%  max b'y - 0.5x'Qx -0.5(x-xo)'RD(x-xo) subject to A'y+s-Qx=c, y free, s>=0,            (2)
%  
%  In each iteration Partial_cholesky is called
%  for solving the linear system of equations GR Dy =(A((Q_tilde)^{-1})A^{T}+RD)=rhs.
%
%  [x,y,s, optimal_solution] = ipdriver(c,A,b,Q,k,maxit_ipm,maxit_Kryl,printlevel,gamma,delta)
%  finds the regularized primal-dual solution x, (y,s) of
%  the programs in form (1)-(2). y are the Lagrange multipliers to the
%  equality constraints and s>=0 the Lagrange multipliers to the
%  nonnegativity constraints on x. optimal_solution indicates if the
%  algorithms solves the problem.
%
%  
%  Dual residual and the duality measure are all below tol.
%  Default: 1e-7.
%
%  maxit: specifies the maximum number of
%  iterations. Default: 100.
%
%  0: turn off iteration output
%  1: print primal and dual residual and duality measure
%  2: print centering parameter and step length
%  3: print residuals in the solution of the step equations
%  Default: 1
%
%  small part of the code was based from the code of Lukas Schork the assignment 1 
%  from the course Large Scale Optimization with Data Science
%
%  Author: Angelis Tzouchas

  [m n] = size(A);
    
  % Make sure that b and c are column vectors of dimension m and n.
  if (size(b,2) > 1) b = b'; end
  if (size(c,2) > 1) c = c'; end
  if (~isequal(size(c),[n,1]) || ~isequal(size(b),[m,1]))
    error('problem dimension incorrect');
  end

  % Make sure that A is sparse and b, c are full.
  if (~issparse(A)) A = sparse(A); end
  if (issparse(b)) b = full(b); end
  if (issparse(c)) c = full(c); end

  % Set default values for missing parameters.
  %if (nargin < 4 || isempty(tol)) tol = 1e-7; end
  if (nargin < 5 || isempty(maxit_ipm)) maxit_ipm = 100; end
  if (nargin < 6 || isempty(printlevel)) printlevel = 1; end
  pl = printlevel;
  
  ep = 1e-2;
  ed = 1e-2;
  eo = 1e-3;
  ekm = 1e-8;
  
  
  Rd = (delta*delta).*speye(m,m);
  
  Rp = (gamma*gamma).*speye(n,n);
  
  Rd_value = delta*delta;
  Rp_value = gamma*gamma;
  
  
  x = ones(n,1);
  y = zeros(m,1);
  s = ones(n,1);  
  
  primal_infeasibility = b-A*x; 
  dual_infeasibility = c-transpose(A)*y-s;
  
  % Set fixed parameters, see below.
  sigma1   = 0.5;
  sigmamin = 0.05;
  sigmamax = 0.95;

  % Initialize persistent variables.
  iter = 0;
  alpha_x = 0;
  alpha_s = 0;
  
  % auxilary vector e
  e = ones(n,1);

  header(pl);
  % initialization of mu = 1  
  mu = transpose(x)*s/n;
  
  optimal_solution = 0;
  
  while (iter < maxit_ipm )

    %from here i took the primal /dual infeasibilities
    primal_infeasibility = b-A*x;
    dual_infeasibility = c-A'*y-s;
    mu_vector = mu*ones(n,1);
    ksmu = mu.*(e) -(x.*s);
    %%%f_double_dash = dual_infeasibility - (1./x).*(mu_vector)+s;
    f_double_dash = dual_infeasibility - (1./x).*ksmu;
    h_double_dash = primal_infeasibility;
    
    %I should check the while condition
    if ((norm(primal_infeasibility,inf )/(1+norm(b, inf)))< ep && (norm(dual_infeasibility, inf)/(1+norm(c, inf))) < ed && (((x')*s/n)/(1+norm((c')*x+(1/2)*(x')*Q*x))) < eo)
        fprintf('optimal solution found\n');
        optimal_solution = 1;
      break;
    end
    
    %if ((norm(primal_infeasibility)/(1+norm(b)))< ep && (norm(dual_infeasibility)/(1+norm(c))) < ed && (((x')*s/n)/(1+norm((c')*x+(1/2)*(x')*Q*x))) < eo)
        %fprintf('optimal solution found\n');
        %optimal_solution = 1;
      %break;
    %end
    
    
    
    iter = iter + 1;
    
    if (iter == 1)
      sigma = sigma1;
    else
      sigma = max(1-alpha_x,1-alpha_s)^5;
    end
    sigma = min(sigma,sigmamax);
    sigma = max(sigma,sigmamin);
    q = sigma*mu*ones(n,1) - x.*s;   
    
    %primal_infeasibility = b-A*x;
    %dual_infeasibility = c-A'*y-s;
    %mu_vector = mu*ones(n,1);
    %ksmu = mu.*(e) -(x.*s);
    %f_double_dash = dual_infeasibility - (1./x).*(mu_vector)+s;
    %f_double_dash = dual_infeasibility - (1./x).*q;
    %h_double_dash = primal_infeasibility;
    
    
    
    
    % update the variables
    theta = x.*(1./s);
    inv_theta = 1./theta;
    Q_tilde = inv_theta + Rp_value;
    inv_Q_tilde = (1./Q_tilde);         
    
    rhs_of_pcg = h_double_dash + A*(inv_Q_tilde.*f_double_dash);
    
    GR_operator = @(x) (A*(inv_Q_tilde.*(A'*x))+Rd*x);
    
    d = sparse(m,1);
    
    for i=1:n
        d = d + (A(:,i).^2)*inv_Q_tilde(i);
    end
    
    d = d + Rd_value;
    
    PS = partial_cholesky_decomposition(k,@(x) GR_operator(x),d);    
    
    [dy, flag,rel_res,Kryl_it] = pcg(@(x) GR_operator(x),rhs_of_pcg,ekm ,maxit_kryl,@(x) inverse_preconditioner(x,PS));
    flag
    rel_res
    Kryl_it
    
    
    
    dx = inv_Q_tilde.*(A'*dy-f_double_dash);
    
    % Q_tilde dx = A'*dy -f_double_dash 
    %ds = (1./x).*(q - s.*dx);
    ds = (1./x).*q - inv_theta.*dx;
    
    
    idx = dx < 0;
    ids = ds < 0;
    alphamax_x = min([1;-x(idx)./dx(idx)]);
    alphamax_s = min([1;-s(ids)./ds(ids)]);
    rho = 0.95;
    alpha_x = rho*alphamax_x;
    alpha_s = rho*alphamax_s;

    % Make step.
    x = x+alpha_x*dx; y = y+alpha_s*dy; s = s+alpha_s*ds;

    % Print iteration output.
    
    xinf = norm(A*x-b);
    sinf = norm(A'*y+s-c);
    %mu = transpose(x)*s/n;
    mu = 0.1*mu;
    output(pl,iter,xinf,sinf,mu,sigma,alpha_x,alpha_s);
    %output(pl,iter,xinf,sinf,mu,sigma,alpha_x,alpha_s,res1,res2,res3);
    
    

  end % while (iter < maxiter)

  % The IPM has terminated either because the solution accuracy is
  % reached or the maximum number of iterations is exceeded. Print
  % result.
  
  fprintf('iterations: %4d\n', iter);
  fprintf('primal feasibility: %8.2e\n', norm(A*x-b));
  fprintf('dual feasibility: %8.2e\n', norm(A'*y+s-c));
  fprintf('complementarity: %8.2e\n', dot(x,s)/n);
end


% ======================================================================
% header
% ======================================================================

function header(pl)
  if (pl >= 1)
    fprintf(' ');
    fprintf('%4s    ', 'iter');
    fprintf('%8s  ', 'pr feas');
    fprintf('%8s  ', 'dl feas');
    fprintf('%8s  ', 'mu');
  end
  if (pl >= 2)
    fprintf('  ');
    fprintf('%8s  ', 'sigma');
    fprintf('%8s  ', 'alpha_x');
    fprintf('%8s  ', 'alpha_s');
  end
  if (pl >= 3)
    fprintf('  ');
    fprintf('%8s  ', 'res1');
    fprintf('%8s  ', 'res2');
    fprintf('%8s  ', 'res3');
  end
  if (pl >= 1)
    fprintf('\n ====    ========  ========  ========');
  end
  if (pl >= 2)
    fprintf('    ========  ========  ========');
  end
  if (pl >= 3)
    fprintf('    ========  ========  ========');
  end
  if (pl >= 1) fprintf('\n'); end
end


% ======================================================================
% output
% ======================================================================

function output(pl,it,xinf,sinf,mu,sigma,alpha_x,alpha_s,res1,res2,res3)
  if (pl >= 1)
    fprintf(' ');
    fprintf('%4d    ', it);
    fprintf('%8.2e  ', xinf);
    fprintf('%8.2e  ', sinf);
    fprintf('%8.2e  ', mu);
  end
  if (pl >= 2)
    fprintf('  ');
    fprintf('%8.2e  ', sigma);
    fprintf('%8.2e  ', alpha_x);
    fprintf('%8.2e  ', alpha_s);
  end
  if (pl >= 3)
    fprintf('  ');
    fprintf('%8.2e  ', res1);
    fprintf('%8.2e  ', res2);
    fprintf('%8.2e  ', res3);
  end
  if (pl >= 1) fprintf('\n'); end
end