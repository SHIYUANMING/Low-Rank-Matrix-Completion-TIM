%function [Xcg, xcost]=RiemOpt_fixedrank(K,r, P, X0)
function [Nres, Xcg, xcost, info, options]=EmbG_fixedrankTR(K,r, P, X0, options)
% Input:
% KxK of rank r matrix
% P: mask matrix
% X0: initial point
%samples: number of samples for all the off diagnoal elements

%Output
%xcost: cost value
%Nres: normalized residual, i.e., ||P(XY-M)||/||P(M)||

    %% Pick the manifold of matrices of size mxn of fixed rank k.
       problem.M = fixedrankembeddedfactory(K, K, r);

    %% Define the problem cost function f(X) = 1/2 * || P.*(X-A) ||^2
    problem.cost = @cost;  % The input X is a structure with fields U, S, V representing a rank k matrix as U*S*V'.
    function f = cost(X)
        Xmat = X.U*X.S*X.V';
        f = .5*norm( P.*Xmat - eye(K), 'fro')^2;
    end

    %% Define the Euclidean gradient of the cost function nabla f(X) = P.*(X-A)
    problem.egrad = @egrad;
    function G = egrad(X)
        % Same comment here about Xmat.
        Xmat = X.U*X.S*X.V';
        G = P.*Xmat - eye(K);
    end


    %% Define the Euclidean Hess  
    problem.ehess = @euclidean_hessian;
    function ehess = euclidean_hessian(X, H)
        ambient_H = problem.M.tangent2ambient(X, H);
        Xdot = ambient_H.U*ambient_H.S*ambient_H.V';
        ehess = P.*Xdot;
    end
    
      %checkhessian(problem);

    %% Notice that for this solver, the Hessian is not needed.
       [Xcg, xcost, info, options] = trustregions(problem, X0, options);
       Nres=sqrt(xcost*2)/sqrt(K);
       
end
