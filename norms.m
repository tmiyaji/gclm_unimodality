function [Nu,z2,C,Ck] = norms(N,R,k,alp,ord,u)
% This MATLAB file computes some bounds for obtained approximate solution.
% Nu(i): encloses (r_{2N}, \phi_i)_{L^2}
% z2: bounds \|(I-P_N)r_{2N}\|_{s-1} from above

  one = intval('1');  two = intval('2'); Pi = intval('pi');

  f = function_f_intval(N,R,k,alp,u);

  v = (intval(1:N).^2)' .* u; % -(u_N)''
  x = intval(zeros(2*N,1));
  x(1:N) = v;
  d = -x + f; % residual (u_N)'' + f(u_N)
  t1 = sup(sqrt(Pi/two)*norm(d));

  Nu = intval(zeros(N,1));
  Nu(1:N) = Pi/two*d(1:N);
  t2 = max(mag(Nu));

  x = d(N+1:2*N);

  z2 = intval(sup(sqrt(Pi/two)*norm(x.*(intval((N+1:2*N)).^(ord-1)))));

  formatSpec = 'defect:%10.4e  |Nu|:%10.4e  z2:%10.4e\n';
  fprintf(formatSpec,t1,t2,sup(z2))

  C = one/intval(N+1);
  Ck = @(k) sqrt(two / ((2*ord-2*k-1)*Pi*intval(N)^(2*ord-2*k-1)));
