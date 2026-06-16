function df = function_df_intval(N,R,alp,j,u)
% compute nonlinear differentiation of f[u]*phi_j by interval arithmetic
% u must be an interval vector
  df = intval(zeros(2*N,1));
  two = intval('2');
  Alp = intval(alp); jot = intval(j);

  % % Code 1
  % % partial derivatives of f(1)
  % if 1<=j-1 && j-1 <= N
  %     df(1) = df(1) - u(j-1)*Alp*(one/intval(j) + one/intval(j-1));
  % end
  % if 1<=j+1 && j+1<=N
  %     df(1) = df(1) - u(j+1)*Alp*(one/intval(j+1) + one/intval(j));
  % end
  % % partial derivatives of f(n) for 2<=n<=N-1
  % for n=2:N-1
  %     en = intval(n);
  %     if 1<=n-j && n-j <=N
  %         df(n) = df(n) + u(n-j)*(-two + Alp*(intval(n-j)/intval(j) + intval(j)/intval(n-j)));
  %     end
  %     if 1<=j-n && j-n<=N
  %         df(n) = df(n) -en*u(j-n)*Alp*(one/intval(j) + one/intval(j-n));
  %     end
  %     if 1<=j+n && j+n<=N
  %         df(n) = df(n) -en*u(j+n)*Alp*(one/intval(j+n)+one/intval(j));
  %     end
  % end
  % % partial derivatives of f(n) for N<=n<=2*N
  % for n=N:2*N
  %     if 1<=n-j && n-j <=N
  %         df(n) = df(n) + u(n-j)*(-two + Alp*(intval(n-j)/intval(j) + intval(j)/intval(n-j)));
  %     end
  % end

  % % Code 2
  n_s = 1;
  n_t = j-1;
  if n_s <= n_t
    n_range = n_s:n_t;
    idx = j - n_range';
    df(n_range) = df(n_range) + Alp * u(idx) .* (intval(idx) / jot - jot ./ intval(idx));
  end

  n_s = 1;
  n_t = N-j;
  if n_s <= n_t
      n_range = n_s:n_t;
      idx = j + n_range';
      df(n_range) = df(n_range) + Alp * u(idx) .* (jot ./ intval(idx) - intval(idx) / jot);
  end

  n_s = max(2,1 + j);
  n_t = min(2*N, N + j);
  if n_s <= n_t
      n_range = n_s:n_t;
      idx = n_range' - j;
      df(n_range) = df(n_range) + u(idx) .* (-two + Alp * (intval(idx)/jot + jot ./ intval(idx)));
  end

  df = intval(R)/two*df;

