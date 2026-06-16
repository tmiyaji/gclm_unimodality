% This MATLAB file computes an approximate solution for gCLM eq.:
%  --------------------------------------------------------------------------  
%      -u^(2)=f(u)=:R*(B(Lu, u) + sin(k*x)) in (0,Pi)
%      B(Lu, u) = -alp*v*u' + v'*u, where v'=Hu (Hilbert transform)
%      u=0 on boundary of (0,Pi)
%  --------------------------------------------------------------------------  
% by Fourier-Galarkin method with usual floating point arithmetic.
  tic
% parameter input
  format long;
  input_data
  e = e_ap;     % stopping criterion (k=3, R>=1000)
  maxit = 100;  % maximum iteration number
  % tol = eps^2;  % tolelance of approximate solution
  u = initial_value(N);      % set initial value
  D = 0.5*pi*diag(1:N)^2;      % diagonal matrix D derived from -u^(2)
  D = full(D);

  it = 1; err = realmax;
  formatSpec = 'R = %-6d N = %-6d alp = %-6d k = %-2d\n';
  fprintf(formatSpec,R,N,alp,k)
  fprintf('%5s %15s\r','its','relative error');
% Newton-Raphson iteration start
  while err > e 
    f = function_f(N,R,k,alp,u);
    g = D*u - 0.5*pi*f(1:N);

    G = D;             % Jacobian matrix G
    for j = 1:N
      df = function_df(N,R,alp,j,u);
      G(:,j) = G(:,j) - 0.5*pi*df(1:N);
    end 
   
    x = u - G\g;              % solve linear system
    % cut small values (Maybe it's better to turn it off for small R. However, tol=eps^3 worked for R=1000)
    if tol > 0
      x(abs(x)<tol) = 0;        
    end
    % if norm(x,inf) < eps      % check some differences
    %   error('obtain trivial zero solution numerically')
    %   u = x;
    % end
    error1 = norm(x-u,inf)/norm(u,inf);
    error2 = norm(g,inf)/norm(u,inf);
    err = max(error1,error2);
    formatSpec = '%5d % 10.5e\n';
    fprintf(formatSpec,it,err)

    % prepare next iteration
    if it < maxit
      it = it + 1; u = x;
    else
      error('not converge: exceed maximum iteration count!')
    end 
  end % end of iteration
 
% write approximate solution to the file as binary
   fid = fopen('appb.dat','w');
   fwrite(fid,u,'double');
   status = fclose(fid);

  toc