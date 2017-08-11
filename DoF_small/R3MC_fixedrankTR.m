%function [Xcg, xcost]=RiemOpt_fixedrank(K,r, P, X0)
function [Xcg, xcost, info, options]=R3MC_fixedrankTR(K,r, P, X0, options)
% Input:
% KxK of rank r matrix
% P: mask matrix
% X0: initial point
%samples: number of samples for all the off diagnoal elements

%Output
%xcost: cost value

    %% Pick the manifold of matrices of size mxn of fixed rank k.
       problem.M = fixedrankfactory_3factors_preconditioned(K, K, r);

    %% Define the problem cost function f(X) = 1/2 * || P.*(X-A) ||^2
    problem.cost = @cost;  % The input X is a structure with fields U, S, V representing a rank k matrix as U*S*V'.
    function f = cost(X)
        Xmat = X.L*X.S*X.R';
        f = .5*norm( P.*Xmat - eye(K), 'fro')^2;
    end

    %% Define the Euclidean gradient of the cost function nabla f(X) = P.*(X-A)
    problem.grad = @(X) problem.M.egrad2rgrad(X, egrad(X));
    function g = egrad(X)
        % Same comment here about Xmat.
        Xmat = X.L*X.S*X.R';
        G = P.*Xmat - eye(K);
        g.L= G*X.R*X.S';
        g.S = X.L'*G*X.R;
        g.R = G'*X.L*X.S;
    end


    %% Define the Euclidean Hess  
    problem.hess = @(X, L) problem.M.ehess2rhess(X, egrad(X), ehess(X, L), L);
    function Hess = ehess(X, eta)
        Xmat = X.L*X.S*X.R';
        G = P.*Xmat - eye(K);
        Pdot  = P.*(eta.L*X.S*X.R' + X.L*eta.S*X.R' + X.L*X.S*eta.R');
        
        Hess.L = G*(X.R*eta.S' + eta.R*X.S') + Pdot*X.R*X.S';
        Hess.R = G'*(X.L*eta.S + eta.L*X.S) + Pdot'*X.L*X.S;
        Hess.S = X.L'*Pdot*X.R + eta.L'*G*X.R + X.L'*G*eta.R;
        
    end
    
%      checkhessian(problem);

    %% Notice that for this solver, the Hessian is not needed.
       [Xcg, xcost, info, options] = trustregions(problem, X0, options);
       
end
