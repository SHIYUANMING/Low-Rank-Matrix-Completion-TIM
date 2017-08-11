function rank=LMaFit_A(K,P)
% K users
% P: # mask for observed entries

%% setting for LMaFit
           opts.tol = 1e-7;
           opts.maxit = 100000;
        
        opts.Zfull = 0;
        opts.est_rank=2; % the increasing rank strategy
        opts.rank_max=K;
        opts.rk_inc=1;

%% Compute an initial guess. Points on the manifold are represented as  structures with three fields: U, S and V. U and V need to be  orthonormal, S needs to be diagonal.
Idx=find(P==1); Ms=eye(K); b=Ms(Idx); % input data for sampling

[X_LMaFit,Y_LMaFit,Out_LMaFit] = lmafit_mc_adp(K,K,1,Idx,b,opts); % estimated rank=1

rank=Out_LMaFit.rank;