
cd ..
cd Dataset;
maxit_ipm =250;
maxit_kryl =400;
Q=0;

gamma = 1e-4;
delta = 1e-3;

Rd_max = 1e-4;
Rd_min = 1e-6;

a = dir;
[p , l] = size(a);

k = 50;

vector_with_solved_problems = zeros(p,1);

vector_which_indicates_free_variable_models = zeros(p-2,1);

for i  = 3: p
    NS = load(a(i).name);
    C = NS.model.obj;
    A = NS.model.A;
    [m,n] = size(A);
    
    if k> m
        k=m;
    end
    
    
    b = NS.model.rhs;
    size(NS.model.obj)


    NS = load(a(i).name)

    c = NS.model.obj;
    A = NS.model.A;
    b = NS.model.rhs;
    sense = NS.model.sense;
    lb = NS.model.lb;
    ub = NS.model.ub;

    cd ..;
    cd matrix-free-ipm;

    [ c, A, b, free_variables, objective_const_term ] = LP_Convert_to_Standard_Form( c, A, b, lb, ub, sense )


    if (isempty(free_variables))
        fprintf('yes')
        vector_which_indicates_free_variable_models(i-2)=1;
        Q = 0;
        %[x,y,s, optimal_solution] = ipm_matrix_free_driver(c,A,b,Q,k,maxit_ipm,maxit_kryl,1,gamma,delta)
        %[x,y,s, optimal_solution] = ipm_matrix_free_driver_dynamic_reg(c,A,b,Q,k,maxit_ipm,maxit_kryl,1,Rd_max,Rd_min,gamma)
        [x,y,s, optimal_solution] = ipm_matrix_free_driver_st(c,A,b,Q,k,maxit_ipm,maxit_kryl,1,gamma,delta)
        vector_with_solved_problems(i-2) = optimal_solution;

    else
        fprintf('no')
        
    end
    cd ..;
    cd Dataset;
    
end