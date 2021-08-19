function [NS]=partial_cholesky_decomposition(k,GR,diagonal)
%  [NS] = partial_cholesky_decomposition(k, GR, diagonal).
%  Function for computing the partial Cholesky decomposition, with k pivots, of 
%  matrix GR with the largest diagonal pivoting.
% 
%  
%  k: indicates the number of pivots for the partial decomposition.
%  
%  GR: indicates the GR matrix where we take advantage of the anonymous
%  functions for performing.
%  
%  diagonal: indicates the diagonal of the GR matrix. 
%
%  Suppose the GR matrix has size m,m
% 
%  NS.L11: k*k      sparse matrix,      
%  NS.L21: (m-k)*k  sparse matrix,    
%  NS.DL vector with the k largest pivots, 
%  NS.DS vector with the diagonal elements of the Schur complement,
%  with dimension m-k.
%
%
    if isa(GR,'function_handle')
        type=0;
    else
        type=1;
        if (~issparse(GR))                                                                                
            GR = sparse(GR);
        end
    end

    if (nargin<3 || isempty(diagonal))
        if type
            diagonal = spdiags(GR,0);  
        else
            error('less number of arguments in the Partial Cholesky decomposition')
        end
    end

    if k>length(diagonal)
        k=length(diagonal);        
    end
    n=length(diagonal);
    I=[];                                                                                                
    Ic=(1:n)';                                                                                           
    P=sparse(k,1);                                                                                       
    B=sparse(n,k);                                                                                       
    NS = struct();                                                                                       
    
    % Computation of the partial Cholesky decomposition
    
    for j=1:k

        [pivot,i] = max(diagonal);                                                                                    
        I=[I;i];
        Ic(i)=0;                                                                                 
        e_i_vector = zeros(n,1);
        e_i_vector(i) = 1;
        
        if type                                                                                          
            b = GR(:,i);
        else
            b = GR(e_i_vector);
        end

        if j>1                                                                                            
            for q=1:j-1
                bq = B(:,q);
                b = b-((1/P(q))*bq(i)).*bq; 
            end
        end
        b(I)=0;
        % Update the diagonal vector.
        diagonal = diagonal - (1/pivot).*(b.^2);          

        diagonal(i) = 0;
        % Store the pivot. 
        P(j) = pivot;
        % Store the Schur complement column. 
        B(:,j) = b;                                                                                        

    end
    
    % Calculation and store of all the parts of the partial Cholesky
    % decomposition.
    
    Ic(Ic==0) = [];
    II = [I;Ic];                                                                                          
    L = B(II,:)./P';                                                                                             
    L_T = L';                                                                                            
    diagonal = diagonal(Ic);
    NS.L11 = (L_T(:,1:k))';
    % Store all the parts of the partial Cholesky decomposition.
    NS.L11 = NS.L11+speye(k);
    NS.L21 = (L_T(:,k+1:end))';              
    NS.DL = P;
    NS.DS = diagonal;
    % Store the permutation vector.
    NS.perm = II;
    
end