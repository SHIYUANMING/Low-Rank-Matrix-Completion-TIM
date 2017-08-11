clc; clear all;

    % Options for the ManOpt (not mandatory)
       options.maxiter = 500;
%     options.maxinner = 30;
%     options.maxtime = 120;
    options.tolgradnorm = 1e-6;
%     options.Delta_bar = min(m, n)*r;
%     options.Delta0 =  min(m, n)*r/64;

%% Problem data
error_th=10^(-6);  % estimation criteria
K=20;  % K-user interference channel
max_realization=100; % # realizations of the network topology
Interfer_link_set = [0,5,10,20,30,40]; 
%Interfer_link_set=[0:10:40];   

R3MC=true; R3MC_SVD=false; EmbG=true; 
R3MCCGRP=true; LMaFit=true;

for loop_interfer=1:length(Interfer_link_set)
    Interfer_link=Interfer_link_set(loop_interfer); % nunber of interference links
    
    rank_temp_EmbG=0;  rank_temp_R3MC=0; rank_temp_R3MC_SVD=0; 
    rank_temp_R3MC_CGRP=0; rank_temp_LMaFit=0;
    
    tsolve_temp_EmbG=0;  tsolve_temp_R3MC=0; tsolve_temp_R3MC_SVD=0; 
    tsolve_temp_R3MC_CGRP=0; tsolve_temp_LMaFit=0;
    
    for loop_realization=1:max_realization
        
    %% Generate a random mask for observed entries: P(i, j) = 1 if the entry  (i, j) of A is observed, and 0 otherwise.
        P = make_rand_Omega(K, Interfer_link);
     
        if R3MC==true
            tstart = clock;
            rank_temp_R3MC=rank_temp_R3MC+R3MC_TRRP(K,P, error_th,options);
            tsolve_temp_R3MC =tsolve_temp_R3MC+etime(clock,tstart);
        end
        
        if R3MC_SVD==true
            tstart = clock;
            rank_temp_R3MC_SVD=rank_temp_R3MC_SVD+R3MC_CG_SVDupdate(K,P, error_th, options); 
            tsolve_temp_R3MC_SVD =tsolve_temp_R3MC_SVD+etime(clock,tstart);
        end
        
        if R3MCCGRP==true
            tstart = clock;
            rank_temp_R3MC_CGRP=rank_temp_R3MC_CGRP+R3MC_CGRP(K,P, error_th, options); 
            tsolve_temp_R3MC_CGRP =tsolve_temp_R3MC_CGRP+etime(clock,tstart);
        end
        
        if EmbG==true
            tstart = clock;
            rank_temp_EmbG=rank_temp_EmbG+EmbG_CGRP(K,P, error_th, options);
            tsolve_temp_EmbG =tsolve_temp_EmbG+etime(clock,tstart);
        end
        
        if LMaFit==true
            tstart = clock;
            rank_temp_LMaFit=rank_temp_LMaFit+LMaFit_A(K,P);
            tsolve_temp_LMaFit =tsolve_temp_LMaFit+etime(clock,tstart);
        end
        
    end
    if R3MC==true
    DoF_R3MC(loop_interfer)=1/(rank_temp_R3MC/max_realization);
    Time_R3MC(loop_interfer)=tsolve_temp_R3MC/max_realization;
    end
    if R3MC_SVD==true
    DoF_R3MC_SVD(loop_interfer)=1/(rank_temp_R3MC_SVD/max_realization);
    Time_R3MC_SVD(loop_interfer)=tsolve_temp_R3MC_SVD/max_realization;
    end
    if R3MCCGRP==true
    DoF_R3MC_CGRP(loop_interfer)=1/(rank_temp_R3MC_CGRP/max_realization);
    Time_R3MC_CGRP(loop_interfer)=tsolve_temp_R3MC_CGRP/max_realization;
    end
    if EmbG==true
    DoF_EmbG(loop_interfer)=1/(rank_temp_EmbG/max_realization);
    Time_EmbG(loop_interfer)=tsolve_temp_EmbG/max_realization;
    end
    if LMaFit==true
    DoF_LMaFit(loop_interfer)=1/(rank_temp_LMaFit/max_realization);
    Time_LMaFit(loop_interfer)=tsolve_temp_LMaFit/max_realization;
    end
end

% DoF_R3MC(1)=1; DoF_R3MC_SVD(1)=1; DoF_EmbG(1)=1;  DoF_R3MC_CGRP(1)=1; DoF_LMaFit(1)=1;

if R3MC==true
plot(Interfer_link_set,DoF_R3MC,'r-o','LineWidth',1.5, 'MarkerSize',10); %R3MC
hold on;
end
if R3MC_SVD==true
plot(Interfer_link_set,DoF_R3MC_SVD,'g-*','LineWidth',1.5, 'MarkerSize',10); %R3MC_SVD
hold on;
end
if R3MCCGRP==true
plot(Interfer_link_set,DoF_R3MC_CGRP,'m-^','LineWidth',1.5, 'MarkerSize',10); %R3MC_SVD
hold on;
end
if EmbG==true
plot(Interfer_link_set,DoF_EmbG,'b-d','LineWidth',1.5, 'MarkerSize',10); %EmbG
hold on;
end
if LMaFit==true
plot(Interfer_link_set,DoF_LMaFit,'k-s','LineWidth',1.5, 'MarkerSize',10); %EmbG
hold on;
end

% save('DoF_R3MC_SVD.mat','DoF_R3MC_SVD');
save('DoF_R3MC.mat','DoF_R3MC'); save('DoF_R3MC_CGRP.mat','DoF_R3MC_CGRP'); 
save('DoF_EmbG.mat','DoF_EmbG'); save('DoF_LMaFit.mat','DoF_LMaFit');

h=legend('R3MC', 'R3MC+CGRP', 'EmbG', 'LMaFit');  
xlabel('Interfering Links Per User','fontsize',12,'fontweight','b','fontname','helvetica');
ylabel('DoF Per User','fontsize',14,'fontweight','b','fontname','helvetica');