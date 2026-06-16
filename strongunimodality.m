%%% Verifiy strong unimodality
tic
intvalinit('displayinfsup');
intvalinit('sharpivmult')  % use sharp interval matrix multiplication
format long;
input_data

if ord < 3
  error('ord must be greater than or equal to 3');
  
end

M = M_su;  % range partition number
Pi = intval('pi');   unit = intval('[0,1]');  zero = intval(zeros(M,1));
two = intval('2');
load('u.dat','-mat');
load('w.dat','-mat'); w=d;
load('beta.dat','-mat'); beta=gamma;

Ck = @(k) sqrt(two / ((2*ord-2*k-1)*Pi*intval(N)^(2*ord-2*k-1)));

I = Pi * (unit + (0:M-1)') / M;

% compute range of u'(x) on each subintervals 
v = u + w;  v = (1:N)' .* v;
k = (1:N)';
R1 = intval(zeros(M, 1));
for i = 1:M
    R1(i) = cos(I(i) * k') * v;
end
R1i = zero;  R1i(:) = intval('[-1,1]')*Ck(1)*beta;
R1 = R1 + R1i;

disp("[1/3]: sign(u'[0,h]) != sign(u'[π-h,π])?");
if R1(1)*R1(M) >= 0
  disp("Failed");
  return;
else
  disp("Passed");
end

disp("[2/3]: {i | 0 \in [u'([x_i])]} is consecutive?");
t1 = find(in(0, R1));
is_consecutive = all(diff(t1) == 1);
if is_consecutive == 0
    disp("Failed");
    return;
else
    disp("Passed");
end

% compute range of u''(x) on each subintervals 
v = u + w;  v = ((1:N)'.^2) .* v;
R2 = intval(zeros(M, 1));
for i = 1:M
    R2(i) = -sin(I(i) * k') * v;
end
R2i = zero; R2i(:) = intval('[-1,1]')*Ck(2)*beta;
R2 = R2 + R2i;

disp("[3/3]: u'' does not change its sign near u'(x)=0?");
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

disp('Strong unimodality is verified!'); 

toc
