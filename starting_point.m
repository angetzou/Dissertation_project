function [x0,y0,s0] = starting_point(c,A,b)
rweight = 0.1;
delta=1e-3;
[m, n] = size(A);
d = zeros(m,1);

for i = 1 : m
    d(i) = sum(A(i,:).^2);
end

d = d + delta;

%inv_prec = 1./d;

A_tr = A';

%coefficient_matrix = A*A_tr + delta*speye(m,m);

ones_vec = ones(n,1);

coefficient_matrix = @(x) (A*(ones_vec.*(A_tr*x))+delta*speye(m,m)*x);

inv_prec = @(x) (1./d).*x;

%x = pcg(coefficient_matrix,b,1e-3,400, diag(inv_prec));

x = pcg(coefficient_matrix,b,1e-3,400, inv_prec);


x = A_tr*x;

rhs_y = A*(c+rweight*x);

%y = pcg(coefficient_matrix,rhs_y,1e-3,400,diag(inv_prec));

y = pcg(coefficient_matrix,rhs_y,1e-3,400,inv_prec);


e = ones(n,1);

% since in our initial value for mu = 1 we derive the s variables in the
% way s = mu*x^{-1}

s = 1./x;

delta_x = max(-1.5*min(x),0);
delta_s = max(-1.5*min(s),0);

product = 0.5*(x+delta_x*e)'*(s+delta_s*e);
delta_x_c = delta_x+product/(sum(s)+n*delta_s);
delta_s_c = delta_s+product/(sum(x)+n*delta_x);

% output
x0 = x + delta_x_c*e;
s0 = s + delta_s_c*e;
y0 = y;

end