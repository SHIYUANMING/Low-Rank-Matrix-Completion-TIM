function rank=EmbG_CGRP(K,P, error_th,options)
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
        [U, S, V] = svds(eye(K), t); X0.U = U; X0.S = S; X0.V = V;
    else


%% Rank-r update
%         X_est=Xcg.U*Xcg.S*Xcg.V';  
%         EGradient=P.*X_est-eye(K); % E Gradient
%         [Ur, Sr, Vr]=svds(-EGradient, 1); 
% 
%         X_guess=X_est+(Ur*Sr*Vr');  % alpha is from wolf search
%         [U_ret, S_ret, V_ret] = svds(X_guess, t); %retraction to the manifold
%         X0.U = U_ret; X0.S = S_ret; X0.V = V_ret;

%% Rank-r update: Riemannian
        X_est=Xcg.U*Xcg.S*Xcg.V';
        [U, S, V]=svd(X_est);   
        EGradient=P.*X_est-eye(K); % E Gradient
        
        PU=U*U'; PU1=eye(K)-U*U'; PV=V*V'; PV1=eye(K)-V*V';
        Sigma1=PU*(-EGradient)*PV+PU1*(-EGradient)*PV+PU*(-EGradient)*PV1;
        [Ur, Sr, Vr]=svds(-EGradient-Sigma1, 1); 
        
        alpha=0.1;
        Sigma2=Ur*Sr*Vr'; 
        Sigma=Sigma1+Sigma2; % search direction
        X_guess=X_est+alpha*Sigma;
        
        [U_ret, S_ret, V_ret] = svds(X_guess+10^(-3)*eye(K), t); %retraction to the manifold
        X0.U = U_ret; X0.S = S_ret; X0.V = V_ret;
           
    end
    
    [Xcg, error_temp]=EmbG_fixedrankCG(K,t, P, X0,options);
    error=sqrt(2*error_temp/K);
    error_record=[error_record, error];
    
    end
    
end

rank=t;