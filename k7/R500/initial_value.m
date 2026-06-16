function u = initial_value(N)
% generate initial guess for approximate solutions
  u = zeros(N,1);
% For k=2 and R=0.4, 0.1*sin(kx) gives a good approximation 
% u(2) = 0.1;
  % u(1) = -sqrt(2.0/3.0); 
  % u(1) = 0.18; 

% Read initial guess from text file
% Just min(N, length(temp)) items are copied

% % Basic solutions
% fid = fopen("data\auto\sol_k2_basic_R1e2.dat", "r"); % N=512
% fid = fopen("data\auto\sol_k2_basic_R3e3.dat", "r");

% % Nontrivial solutions
% fid = fopen("data\auto\sol_k1_R1e3.dat", "r");
% fid = fopen("data\auto\sol_k2_R1e1.dat", "r"); % N=128
% fid = fopen("data\auto\sol_k2_R1e3.dat", "r"); % R=100,N=256
% fid = fopen("data\auto\sol_k2_R3e3.dat", "r");
% fid = fopen("data\auto\sol_k3_R3e3.dat", "r");
% fid = fopen("data\auto\sol_k4_R5e3.dat", "r");
% fid = fopen("data\auto\sol_k5_R5e3.dat", "r");
% fid = fopen("data\auto\sol_k6_R5e3.dat", "r");
fid = fopen("data\auto\sol_k7_R5e3.dat", "r");
% 
data = textscan(fid, '%f');
fclose(fid);
temp = data{1}(2:end);
len_temp = length(temp);
u(1:min(N,len_temp)) = temp(1:min(N,len_temp));
