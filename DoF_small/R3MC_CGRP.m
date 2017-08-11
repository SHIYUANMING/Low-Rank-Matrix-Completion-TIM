function rank=R3MC_CGRP(K,P, error_th,options)
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
%         X_est=Xcg.L*Xcg.S*Xcg.R';
%         [U1, S1, V1]=svds(P.*X_est-eye(K), 1);  X_guess=X_est-U1*S1*V1';
%         [U, S, V] = svds(X_guess+10^(-3)*eye(K), t); X0.L = U; X0.S = S; X0.R = V;

%% Rank-r update: Riemannian
        X.L=Xcg.L; X.S=Xcg.S; X.R=Xcg.R;
        
        Xmat = X.L*X.S*X.R';  % Euclidean Gradient
        G = P.*Xmat - eye(K);
        g.L= G*X.R*X.S';
        g.S = X.L'*G*X.R;
        g.R = G'*X.L*X.S;
        
        rgrad=egrad2rgrad(X, g);  % Riemannian Gradient in vectors
        Sigma1=rgrad.L*X.S*X.R'+X.L*rgrad.S*X.R'+X.L*X.S*rgrad.R';
        
        [Ur, Sr, Vr]=svds(-G+Sigma1, 1); 
        
        alpha=0.1;
        Sigma2=Ur*Sr*Vr'; 
        Sigma=Sigma2-Sigma1; % search direction
        X_guess=Xmat+alpha*Sigma;
        
        [U_ret, S_ret, V_ret] = svds(X_guess+10^(-3)*eye(K), t); %retraction to the manifold
        X0.L = U_ret; X0.S = S_ret; X0.R = V_ret;
           
    end
    
    [Xcg, error_temp]=R3MC_fixedrankCG(K,t, P, X0,options);
    error=sqrt(2*error_temp/K);
    error_record=[error_record, error];
    
    end
    
end

rank=t;

end


