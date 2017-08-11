function rank=R3MC_CG_SVDupdate(K,P, error_th,options)
% K users
% Interfer_link: # interfering links


%% Compute an initial guess. Points on the manifold are represented as  structures with three fields: U, S and V. U and V need to be  orthonormal, S needs to be diagonal.
error=10^(10); t=0; error_record=[];
while error>=error_th
    t=t+1;
    if t==K
        break;
    else
        
    if t==1
        [U, S, V] = svds(eye(K), t); X0.L = U; X0.S = S; X0.R = V;
    else

%% Rank-one update
        X_est=Xcg.L*Xcg.S*Xcg.R';
        [U1, S1, V1]=svds(P.*X_est-eye(K), 1);  X_guess=X_est-U1*S1*V1';
        [U, S, V] = svds(X_guess+10^(-3)*eye(K), t); X0.L = U; X0.S = S; X0.R = V;
           
    end
    
    [Xcg, error_temp]=R3MC_fixedrankCG(K,t, P, X0,options);
    error=sqrt(2*error_temp/K);
    error_record=[error_record, error];
    
    end
    
end

rank=t;