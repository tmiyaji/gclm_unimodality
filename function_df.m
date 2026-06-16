function df = function_df(N,R,alp,j,u)
% compute derivative of f w.r.t. u(j) by floating point arithmetic
% 1 <= j <= N 
  df = zeros(2*N,1);
  n_s = 1;
  n_t = j-1;
  if n_s <= n_t
    n_range = n_s:n_t;
    idx = j - n_range';
    df(n_range) = df(n_range) + alp * u(idx) .* (idx / j - j ./ idx);
  end

  n_s = 1;
  n_t = N-j;
  if n_s <= n_t
      n_range = n_s:n_t;
      idx = j + n_range';
      df(n_range) = df(n_range) + alp * u(idx) .* (j ./ idx - idx / j);
  end

  n_s = max(2,1 + j);
  n_t = min(2*N, N + j);
  if n_s <= n_t
      n_range = n_s:n_t;
      idx = n_range' - j;
      df(n_range) = df(n_range) + u(idx) .* (-2 + alp * (idx/j + j ./ idx));
  end

  df = R/2*df;
