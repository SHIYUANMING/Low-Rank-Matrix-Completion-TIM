function [Omega_sparse] = make_rand_Omega_squarefixed(n,nb_samples)
% Make a random sampling using specified number of samples for the off diagonal entries in nxn matrix.
% Output: Omega_sparse as a sparse matrix (optional output).

Omega_sel=[1:n^2];
Diag_index=n*([1:n]-1)+[1:n];
Omega_sel(Diag_index)=[];  % index for all the off-diagonal elements

if nb_samples==0
    [Omega_i, Omega_j] = ind2sub([n,n], Diag_index);
    Omega_sparse = sparse(Omega_i, Omega_j, 1, n, n, nb_samples+n);
else
    randsampling_index=randperm(n^2-n, nb_samples);
    Omega_sampling=Omega_sel(randsampling_index); % random sampling for the index of all the off-diagonal elements
    
    Omega = [Omega_sampling, Diag_index];  % final index set
    Omega=unique(Omega);
    
    [Omega_i, Omega_j] = ind2sub([n,n], Omega);
    Omega_sparse = sparse(Omega_i, Omega_j, 1, n, n, nb_samples+n);
end


%% Jafar Setting Fig.4, nb_samples=9, n=5
%     Omega_i=[ones(3,1);2*ones(3,1); 3*ones(3,1); 4*ones(3,1); 5*ones(2,1)];
%     Omega_j=[[1;3;4];[2;3;4];[1;3;5];[1;4;5];[2;5]];
%     Omega_sparse = sparse(Omega_i, Omega_j, 1, n, n, nb_samples+n);


%% Jafar Setting Fig.12(a), nb_samples=7, n=4
%     Omega_i=[ones(2,1);2*ones(2,1); 3*ones(3,1); 4*ones(4,1)];
%     Omega_j=[[1;3];[1;2];[1;2;3];[1;2;3;4]];
%     Omega_sparse = sparse(Omega_i, Omega_j, 1, n, n, nb_samples+n);

%% Jafar Setting Fig.12(b), nb_samples=6, n=4
%     Omega_i=[ones(2,1);2*ones(2,1); 3*ones(2,1); 4*ones(4,1)];
%     Omega_j=[[1;3];[1;2];[2;3];[1;2;3;4]];
%     Omega_sparse = sparse(Omega_i, Omega_j, 1, n, n, nb_samples+n);

%% Jafar Setting Fig.10, nb_samples=13, n=7
%     Omega_i=[ones(2,1);2*ones(3,1); 3*ones(3,1); 4*ones(3,1); 5*ones(3,1); 6*ones(3,1); 7*ones(3,1)];
%     Omega_j=[[1;5];[2;4;5];[2;3;6];[3;4;7];[1;2;5]; [1;3;6]; [1;4;7]];
%     Omega_sparse = sparse(Omega_i, Omega_j, 1, n, n, nb_samples+n);

%% Jafar Setting Fig.10, nb_samples=15, n=10
%     Omega_i=[ones(2,1); 2*ones(4,1); 3*ones(3,1); 4*ones(3,1); 5*ones(3,1); 6*ones(2,1); 7*ones(2,1); 8*ones(2,1); 9*ones(3,1); 10*ones(1,1)];
%     Omega_j=[[1;3]; [1;2;5;10]; [2;3;7]; [4;5;10]; [4;5;8]; [3;6]; [3;7]; [1;8]; [6;8;9]; [10]];
%     Omega_sparse = sparse(Omega_i, Omega_j, 1, n, n, nb_samples+n);

end
