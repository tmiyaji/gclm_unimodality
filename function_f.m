function f = function_f(N,R,k,alp,u)
% compute nonlinear function f by floating point arithmetic
  f = zeros(2*N,1);
  p = (u * u') .* (alp * (1:N) ./ (1:N)' - 1);
  q = (u * u') .* (alp * (1:N) ./ (1:N)' + 1);

  % coefficients of sin(m+n)x terms (2..2N)
  xi = zeros(2*N,1);
  for n = 2:2*N
    i_range = max(1, n-N):min(N, n-1);
    j_range = n - i_range;
    idx = sub2ind([N, N], i_range, j_range);
    xi(n) = sum(p(idx));
  end
  % coeffcients of sin(m-n)x terms (-(N-1)..(N-1))
  eta = zeros(2*N-1,1);
  for n = -(N-1):(N-1)
    eta(n + N) = sum(diag(q, -n));
  end
  f(1:N-1) = f(1:N-1) + eta(N+1:2*N-1) - eta(N-1:-1:1);
  f(2:2*N) = f(2:2*N) + xi(2:2*N);

  f = R/2*f; % "/2" comes from the sin addition formula 
  f(k) = f(k) + R;
