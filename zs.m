function [z0,z1,z3,z4,z5,z6,Mu] = zs(N,u,w,alp,C,Ck,ord)
% This MATLAB file computes some bounds in FN-Int

  zero=intval('0'); one= intval('1');  two = intval('2');
  four = intval('4');
  Alp= intval(alp); Pi = intval('pi'); s = sqrt(Pi/two);
  v = u + w;
% -=-=-=-=-= z0 -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= %
% z0(m) encloses \| (I-P_N)(-(a+1)v'_N \phi_m - v_N \phi'_m ) \|_{L^2}
  z0 = intval(zeros(N,1));  
  for m=1:N
    em = intval(m);
    en = intval((1:N)');
    p = ((Alp + one) * en + em) .* v(1:N) / two;
    q = ((Alp + one) * en - em) .* v(1:N) / two;

    kappa = intval(zeros(2*N,1));
    kappa(m+1:N+m) = kappa(m+1:N+m) + p;
    % n from 1-m to -1
    if m > 1
      kappa(1:m-1) = kappa(1:m-1) - q(m-1:-1:1);
    end
    % n from 1 to N-m
    if m < N
      kappa(1:N-m) = kappa(1:N-m) + q(m+1:N);
    end
    kappa(1:N) = zero; % Projection by (I-P_N)
    z0(m) = intval(sup(s*norm(kappa)));
  end

% -=-=-=-=-= z1 -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= %
  z1 = intval(zeros(N,1));    
  for m=1:N
    em = intval(m);
    en = intval((1:N)');
    p = ((Alp + one) * en + two*Alp*em) .* v(1:N) ./ (two*en);
    q = ((Alp + one) * en - two*Alp*em) .* v(1:N) ./ (two*en);

    kappa = intval(zeros(2*N,1));
    % sin(n+m)x coefficients: n from m+1 to m+N
    kappa(m+1:N+m) = kappa(m+1:N+m) + p;
    % sin(n-m)x coefficients
    % n from 1-m to -1
    if m > 1
      kappa(1:m-1) = kappa(1:m-1) - q(m-1:-1:1);
    end
    % n from 1 to N-m
    if m < N
      kappa(1:N-m) = kappa(1:N-m) + q(m+1:N);
    end
    kappa(1:N) = zero; % Projection by (I-P_N)
    z1(m) = intval(sup(s*norm(kappa)));
  end


% -=-=-=-=-= Mu -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= %
% enclosure of (B(Lw_N, w_N),\phi_m)_{L^2}
  Mu = intval(zeros(N,1));
  em = intval((1:N));
  p = (w * w') .* (Alp*em ./ em' - one);
  q = (w * w') .* (Alp*em ./ em' + one);
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
  Mu(1:N-1) = Mu(1:N-1) + eta(N+1:2*N-1) - eta(N-1:-1:1);
  Mu(2:N) = Mu(2:N) + xi(2:N);
  Mu = Mu * Pi / four;

% -=-=-=-=-= z4, z5, z6 -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= %
  en = intval((1:N)');  
  gma = @(k) intval(sup(norm((en.^k).*v, 1)));

  z4_tmp = @(k) nchoosek(ord-1,k) * (abs(Alp) * C^k * gma(k-1) + C^(k+1) * gma(k));
  z4 = zero;
  for k = 0:ord-1
      z4 = z4 + z4_tmp(k);
  end
  z4 = intval(sup(z4));

  z5_tmp = @(k) nchoosek(ord-1,k) * (abs(Alp) * C^(k+2) * gma(k+1) + C^(k+1) * gma(k));
  z5 = zero;
  for k = 0:ord-1
      z5 = z5 + z5_tmp(k);
  end
  z5 = intval(sup(z5));

  z6_tmp = @(k) nchoosek(ord-1,k) * (abs(Alp) * Ck(k-1) * C^k + Ck(k) * C^(k+1));
  z6 = zero;
  for k = 0:ord-1
      z6 = z6 + z6_tmp(k);
  end
  z6 = intval(sup(z6));


% -=-=-=-=-= z3 -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= %
% z3: bounds \|(I-P_N)(B(Lv_N, w_N) + B(Lw_N, u_N))\|_{s-1}
  e = intval(zeros(2*N,1));
  em = intval((1:N));
  p = (v * w' + w * u') .* (Alp*em ./ em' - one);
  q = (v * w' + w * u') .* (Alp*em ./ em' + one);
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
  e(1:N-1) = e(1:N-1) + eta(N+1:2*N-1) - eta(N-1:-1:1);
  e(2:2*N) = e(2:2*N) + xi(2:2*N);
  e = e / two;
  e(1:N) = zero;
  z3 = intval(sup(s*norm(e.*(intval((1:2*N)').^(ord-1)))));
