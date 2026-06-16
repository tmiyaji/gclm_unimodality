%%% Verifiy unimodality

tic
intvalinit('displayinfsup');
intvalinit('sharpivmult')  % use sharp interval matrix multiplication
format long;
input_data

M = M_wu;  % range partition number
Pi = intval('pi');   unit = intval('[0,1]');  zero = intval(zeros(M,1));
two = intval('2');
load('u.dat','-mat');
load('w.dat','-mat'); w=d;
load('beta.dat','-mat'); beta=gamma;

%%% TODO: fix the following
Ck = @(k) sqrt(two / ((2*ord-2*k-1)*Pi*intval(N)^(2*ord-2*k-1)));

% generate M-th subintervals on (0,Pi)
I = zero;
for i = 1: M
I(i) = Pi*(unit + i - 1)/M;
end

% compute range of psi'(x) on each subintervals 
v = u + w;  v = v ./ (1:N)';
k = (1:N)';
R1 = cos(k * I')' * v;
R1i = zero;  R1i(:) = intval('[-1,1]')*Ck(-1)*beta;
R1 = R1 + R1i;

disp("[1/3]: sign(ψ'[0,h]) != sign(ψ'[π-h,π])?");
if R1(1)*R1(M) >= 0
  disp("Failed");
  return;
else
  disp("Passed");
end

disp("[2/3]: {i | 0 \in [ψ'([x_i])]} is consecutive?");
t1 = find(in(0, R1));
is_consecutive = all(diff(t1) == 1);
if is_consecutive == 0
    disp("Failed");
    return;
else
    disp("Passed");
end

% compute range of u''(x) on each subintervals 
v = u + w;
R2 = -sin(k * I')' * v;
R2i = zero; R2i(:) = intval('[-1,1]')*Ck(0)*beta;
R2 = R2 + R2i;


disp("[3/3]: ψ'' does not change its sign near ψ'(x)=0?");
if R1(1) > 0
    if all(R2(t1)<0)
        disp("Passed");
    else
        disp("Failed");
        return;
    end
else
    if all(R2(t1)>0)
        disp("Passed");
    else
        disp("Failed");
        return;
    end
end

disp('Unimodality is verified!'); 


toc
