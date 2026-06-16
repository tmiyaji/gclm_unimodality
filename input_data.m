N = 5000;       % Truncation wave number
R = 1000.0;      % Reynolds number
k = 2;         % Wave number of external force
alp = -2;      % gCLM parameter
ord = 3;       % regularity parameter: X^s (s=ord)
tol = eps^3;   % tolelance of approximate solution (if tol>0, then x(abs(x)<tol) = 0)
e_ap = 10^(-12);% Stopping criterion for the Newton method
M_su = 10000;  % Number of subinterval used in strongunimodality.m
M_wu = 10000;  % Number of subinterval used in unimodality.m
