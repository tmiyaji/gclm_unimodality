function f = function_f_intval(N,R,k,alp,u)
% compute nonlinear function f by interval arithmetic
  Alp = intval(alp); one = intval('1');
  f = intval(zeros(2*N,1));
  em = intval((1:N));
  p = (u * u') .* (Alp * em ./ em' - one);
  q = (u * u') .* (Alp * em ./ em' + one);
  % coefficients of sin(m+n)x terms (2..2N)
  xi = intval(zeros(2*N,1));
  for n = 2:2*N
    i_range = max(1, n-N):min(N, n-1);
    j_range = n - i_range;
    idx = sub2ind([N, N], i_range, j_range);
    xi(n) = sum(p(idx));
  end
  % coeffcients of sin(m-n)x terms (-(N-1)..(N-1))
  eta = intval(zeros(2*N-1,1));
  for n = -(N-1):(N-1)
    eta(n + N) = sum(diag(q, -n));
  end
  f(1:N-1) = f(1:N-1) + eta(N+1:2*N-1) - eta(N-1:-1:1);
  f(2:2*N) = f(2:2*N) + xi(2:2*N);

  f = intval(R)/intval('2')*f; % "/2" comes from the sin addition formula 
  f(k) = f(k) + intval(R);



