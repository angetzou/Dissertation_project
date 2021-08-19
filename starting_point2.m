function [x0,y0,s0] = starting_point2(c,A,b)

rweight = 0.1;
delta=1e-3;

[m, n] = size(A);

d = zeros(m,1);


for i = 1 : m
    d(i) = sum(A(i,:).^2);
end

d = d + delta;

A_tr = A';

ones_vec = ones(n,1);

coefficient_matrix = @(x) (A*(ones_vec.*(A_tr*x))+delta*speye(m,m)*x);

inv_prec = @(x) (1./d).*x;

x = pcg(coefficient_matrix,b,1e-3,400, inv_prec);

x = A_tr*x;

rhs_y = A*(c+rweight*x);

y = pcg(coefficient_matrix,rhs_y,1e-3,400,inv_prec);

e = ones(n,1);

for i = 1:n
    if x(i)<=0
        x(i)=1;
    end
end

% since in our initial value for mu = 1 we derive the s variables in the
% following way
s = 1./x;
s0 = s;
y0 = y;
x0 = x;


end