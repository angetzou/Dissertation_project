function y = inverse_preconditioner(vector,NS)
        
        deperm(NS.perm)=1:length(NS.perm);                                                                              

        n=length(NS.DL);
        vector=vector(NS.perm);
        v1=vector(1:n);
        v2=vector(n+1:end);

        w1=NS.L11\v1;
        w2=v2-NS.L21*w1;

        z1=w1./NS.DL;
        z2=w2./NS.DS;

        y2=z2;
        y1=z1-NS.L21'*y2;
        y1=NS.L11'\y1;

        y=[y1;y2];
        y=y(deperm);
    
end