clc; clear all;

    % Options for the ManOpt (not mandatory)
    options.maxiter = inf;
%     options.maxinner = 30;
%     options.maxtime = 120;
    options.tolgradnorm = 1e-6;
%     options.Delta_bar = min(m, n)*r;
%     options.Delta0 =  min(m, n)*r/64;

R3MC_CG=true; R3MC_TR=false; 
EmbG_CG=false; LMaFit=false;
savedata=false;
%% Problem data
% K=30;  % K-user interference channel
% Interfer_link=floor(4*K); % nunber of interference links

K=100;  % K-user interference channel
Interfer_link=400; % nunber of interference links

%% Generate a random mask for observed entries: P(i, j) = 1 if the entry  (i, j) of A is observed, and 0 otherwise.
P = make_rand_Omega(K, Interfer_link);
% load('P.mat');
Idx=find(P==1); Ms=eye(K); b=Ms(Idx); % input data for sampling

%% Compute an initial guess. 
t=5; [L, S, R] = svds(eye(K), t); 
X0_R3MC.L = L; X0_R3MC.S = S; X0_R3MC.R = R;  % R3MC: based on palor factorization X=L*S*R'
X0_EmbG.U = L; X0_EmbG.S = S; X0_EmbG.V = R;  % EmbG: based on SVD factorization X=U*S*V'

%% setting for LMaFit
        opts.tol = 1e-10;
        opts.maxit = 50000;
        opts.Zfull = 0;
        opts.est_rank=0; % the increasing rank strategy

%% Algorithm: R3MC
if R3MC_CG==true
[Nres_R3MC_CG, Xcg_R3MC_CG, error_temp_R3MC_CG, info_R3MC_CG, option_R3MC_CG]=R3MC_fixedrankCG(K,t, P, X0_R3MC, options);
end

if R3MC_TR==true
[Nres_R3MC_TR, Xcg_R3MC_TR, error_temp_R3MC_TR, info_R3MC_TR, option_R3MC_TR]=R3MC_fixedrankTR(K,t, P, X0_R3MC, options); 
end

%% Algorihtm: EmbG
if EmbG_CG==true
[Nres_EmbG_CG, Xcg_EmbG_CG, error_temp_EmbG_CG, info_EmbG_CG, option_EmbG_CG]=EmbG_fixedrankCG(K,t, P, X0_EmbG, options);
end

%% Algorithm: LMaFit
if LMaFit==true
tstart = clock;
[X_LMaFit,Y_LMaFit,Out_LMaFit] = lmafit_mc_adp(K,K,t,Idx,b,opts);
Time4 = etime(clock,tstart);
end 

figure;
%% R3MC
if R3MC_CG==true
NR1=[info_R3MC_CG.cost];  Time1=[info_R3MC_CG.time]; Time1=Time1(length(Time1));
semilogy([1:length(NR1)]/length(NR1)*Time1, sqrt(NR1.*(2/K)), 'r-','LineWidth',1.2, 'MarkerSize',6);
hold on;
end

if R3MC_TR==true
NR2=[info_R3MC_TR.cost]; Time2=[info_R3MC_TR.time]; Time2=Time2(length(Time2));
semilogy([1:length(NR2)]/length(NR2)*Time2, sqrt(NR2.*(2/K)), 'b-*','LineWidth',1.2, 'MarkerSize',6);
hold on;
end

%% EmbG
if EmbG_CG
NR3=[info_EmbG_CG.cost]; Time3=[info_EmbG_CG.time]; Time3=Time3(length(Time3));
semilogy([1:length(NR3)]/length(NR3)*Time3, sqrt(NR3.*(2/K)), 'g-s','LineWidth',1.2, 'MarkerSize',6);
hold on;
end

%% LMaFit
if LMaFit==true
NR4=Out_LMaFit.obj;
semilogy([1:length(NR4)]/length(NR4)*Time4, NR4, 'm-d','LineWidth',1.2, 'MarkerSize',6);
hold on;
end

if savedata==true
save('/Users/Yuanming/Dropbox/Paper/[8]TIM_LRMC/Simulation/R3MC/Convergence_fixedrank/save/NR.mat','NR1', 'NR2','NR3','NR4','Time1','Time2','Time3','Time4');
end

legend('R3MC with CG', 'R3MC with TR', 'LRGeom', 'LMaFit')
xlabel('Time in Seconds');
ylabel('Normalized Residual');
set(gca,'FontSize',14)