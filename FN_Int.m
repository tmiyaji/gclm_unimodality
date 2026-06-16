% This MATLAB file try to verify the exact solution from obtained
% approximate solution for gCLM eq.:
%  --------------------------------------------------------------------------  
%      -u^(2)=f(u)=:R*(B(Lu, u) + sin(k*x)) in (0,2*Pi)
%      B(Lu, u) = -alp*v*u' + v'*u, where v'=Hu (Hilbert transform)
%      u=0 on boundary of (0,Pi)
%  --------------------------------------------------------------------------  
% by Fourier-Galarkin method with interval arithmetic.
  tic
  intvalinit('displayinfsup');
  intvalinit('sharpivmult')  % use sharp interval matrix multiplication
  format long;
  input_data
  formatSpec = 'R = %-6d N = %-6d alp = %-6d k = %-2d ord = %-2d\n';
  fprintf(formatSpec,R,N,alp,k,ord)

  one = intval('1');  two = intval('2');  unit = intval('[-1,1]');
  Pi = intval('pi'); 
  Alp = intval(alp);

  fid = fopen('appb.dat');
  u = fread(fid,'double');
  status = fclose(fid);
  u = intval(u);

  D = Pi/two*intval(diag(1:N))^2;  % diagonal matrix D
  G = full(D);                 % Jacobian matrix G
  for j = 1:N
    df = function_df_intval(N,R,alp,j,u);
    G(:,j) = G(:,j) - (Pi/two)*df(1:N);
  end 

  [Nu,z2,C,Ck] = norms(N,R,k,alp,ord,u);
  MAXIT = 10;  d0 = 1d-1;
  w = unit*eps*ones(N,1);  beta = eps;     % initial candidate set           
  fprintf('%5s %6s %8s %6s %8s %8s %8s %8s %10s %8s %8s\r',...
	  'its','|z0|','|z1|','z2','z3','z4','z5','z6','|Wn|','|W*|','beta');
  formatSpec = '%5d %7.2e %7.2e %7.2e %7.2e %7.2e %7.2e %7.2e %7.2e %7.2e %7.2e\n';
  for j = 1:MAXIT
    w = (1+d0)*w;  beta = (1+d0)*beta;    % inflation of candidate set
    [z0,z1,z3,z4,z5,z6,Mu] = zs(N,u,w,alp,C,Ck,ord);
    d = R*(unit*(beta*(C^ord)*(z0*C + z1 + beta*(abs(Alp)+one)*C^ord)) + Mu) + Nu;
    d = G\d;                                % next Wn
    gamma = sup(C * (z2 + R * (z3 + beta * (z4 + z5 + z6 * beta))));
    fprintf(formatSpec,...
            j,max(mag(z0)),max(mag(z1)),sup(z2),sup(z3),sup(z4),...
	        sup(z5),sup(z6),max(mag(d)),gamma,beta)
    if gamma > 1
      error('Verification fails!')
    end 
    res = sum(in(d,w));                     % checking contraction: finite
    if gamma <= beta                        % infinite part 
      res = res + 1;
    end 
    if res == N+1                           % verification check
      disp('Verification has been completed!'); 
      fprintf('max(Wn): %d\n',max(mag(d))); % maximum-norm of Wn
      fprintf('  alpha: %d\n',gamma);        % X-norm of W*
      save('u.dat','u')                     % save approximate solution
      save('w.dat','d')                     % save Wn
      save('beta.dat','gamma')              % save W*
      break;
    end
    w = d;  beta = gamma;                   % preparing next step
  end
toc