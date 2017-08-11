%function [Xcg, xcost]=RiemOpt_fixedrank(K,r, P, X0)
function [Nres, Xcg, xcost, info, options]=R3MC_fixedrankCG(K,r, P, X0, options)
% Input:
% KxK of rank r matrix
% P: mask matrix
% X0: initial point
%samples: number of samples for all the off diagnoal elements

%Output
%xcost: cost value
%Nres: normalized residual, i.e., ||P(XY-M)||/||P(M)||

    %% Pick the manifold of matrices of size mxn of fixed rank k.
       problem.M = fixedrankfactory_3factors_preconditioned(K, K, r);

    %% Define the problem cost function f(X) = 1/2 * || P.*(X-A) ||^2
    problem.cost = @cost;  % The input X is a structure with fields L, S, R representing a rank k matrix as L*S*R'.
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
    
      %checkgradient(problem);

    %% Notice that for this solver, the Hessian is not needed.
       [Xcg, xcost, info, options] = conjugategradient(problem, X0, options);
       Nres=sqrt(xcost*2)/sqrt(K);
       
end
