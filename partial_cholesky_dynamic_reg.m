function [NS]=partial_cholesky_dynamic_reg(k,GR,diagonal,Rd_max,Rd_min)
%  [NS] = partial_cholesky_decomposition(k, GR, diagonal,Rd_max,Rd_min).
%  Function for computing the partial Cholesky decomposition, with k pivots, of 
%  matrix GR with the largest diagonal pivoting, for dynamic regularization
% 
%  
%  k: indicates the number of pivots for the partial decomposition.
%  
%  GR: indicates the GR matrix where we take advantage of the anonymous
%  functions for performing.
%  
%  diagonal: indicates the diagonal of the GR matrix. 
%
%  Rd_max: The strong regularization term for pivots which falls dangerously close to zero.
%
%  Rd_min: The minimum regularization term for pivots.
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
    
    % Computation of the Partial Cholesky decomposition with dynamic regularization of the matrix Gr.
    
    indicator_with_type_of_reg = zeros(n,1);
    
    % the respective vector with the type of regularization
    
    for j=1:k

        [pivot,i]=max(diagonal);                                                                                    
        I=[I;i];Ic(i)=0;                                                                                 
        e_i_vector = zeros(n,1);
        e_i_vector(i)=1;
        
        if j == 1
            if type
                b=GR(:,i);
                if pivot < 1e-6
                    pivot = pivot + Rd_max;
                    
                    indicator_with_type_of_reg(i)=1;
                else
                    pivot = pivot + Rd_min;
                    
                end
            else
                b=GR(e_i_vector);
                if pivot < 1e-6
                    pivot = pivot + Rd_max;
                    indicator_with_type_of_reg(i)=1;
                else
                    pivot = pivot + Rd_min;
                    
                end
            end
        end
        
        if j ~= 1
            if type
                b=GR(:,i);
                if pivot < 1e-6
                    pivot = pivot + Rd_max;
                    indicator_with_type_of_reg(i)=1;
                else
                    pivot = pivot + Rd_min;
                end
            else
                b=GR(e_i_vector);
                if pivot < 1e-6
                    pivot = pivot + Rd_max;
                    indicator_with_type_of_reg(i)=1;
                else
                    pivot = pivot + Rd_min;
                end
            end
        end

        if j > 1                                                                                           
            for q=1:j-1
                bq=B(:,q);
                b=b-((1/P(q))*bq(i)).*bq;
            end
        end
        b(I)=0;
        % update the diagonal vector
        diagonal=diagonal-(1/pivot).*(b.^2);          

        diagonal(i)=0;
        % Store the pivot
        P(j)=pivot;
        % Store the schur complement column
        B(:,j)=b;                                                                                        

    end
    % Compute and store all the parts of the partial Cholesky Decomposition the individual blocks of the decomposition LDL^T.
    Ic(Ic==0)=[];
    II=[I;Ic];                                                                                          
    L=B(II,:)./P';                                                                                            
    L_T = L';                                                                                            
    diagonal=diagonal(Ic);
    NS.L11=(L_T(:,1:k))';
    % store all the parts of the Partial Choleksy Decomposition
    NS.L11=NS.L11+speye(k);
    NS.L21=(L_T(:,k+1:end))';              
    NS.DL=P;
    NS.DS=diagonal;
    NS.perm=II;
    % store the permutation vector.
    NS.indicator = indicator_with_type_of_reg;
    %store the vector which indicates the type of regularization
end