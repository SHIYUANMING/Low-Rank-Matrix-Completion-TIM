function rgrad = egrad2rgrad(X, egrad)


    function X = prepare(X)
        if ~all(isfield(X,{'StS','SSt'}) == 1)
            X.SSt = X.S*X.S';
            X.StS = X.S'*X.S;
        end
    end


    symm = @(X) .5*(X+X');

    
        X = prepare(X);
        
        SSL = X.SSt;
        ASL = 2*symm(SSL*(egrad.S*X.S'));
        
        SSR = X.StS;
        ASR = 2*symm(SSR*(egrad.S'*X.S));
        
        [BL, BR] = tangent_space_lyap(X.S, ASL, ASR); % It computes the solution without calling Matlab's Lyap.
        
        rgrad.L = (egrad.L - X.L*BL)/X.SSt;
        rgrad.R = (egrad.R - X.R*BR)/X.StS;
        rgrad.S = egrad.S;
        
        % Debug
        %         BL1 = lyap(SSL, -ASL); % Alternate way
        %         BR1 = lyap(SSR, -ASR);
        %         norm(skew(X.SSt*(rgrad.L'*X.L) + rgrad.S*X.S'), 'fro')
        %         norm(skew(X.StS*(rgrad.R'*X.R) - X.S'*rgrad.S), 'fro')
        
        function[BU, BV] = tangent_space_lyap(R, E, F)
    % We intent to solve a linear system    RR^T  BU + BU RR^T  = E
    %                                       R^T R BV + BV R^T R = F
    % for BU and BV.
    %
    % This can be solved using two calls to the Matlab's lyap.
    % However, we can still have a more efficient implementation
    % that does not require the full functionaliyt of Matlab's lyap.
    
    [U, Sigma, V] = svd(R);
    E_mod = U'*E*U;
    F_mod = V'*F*V;
    b1 = E_mod(:);
    b2 = F_mod(:);
    
    r = size(Sigma, 1);
    sig = diag(Sigma); % all the singular values in a vector
    sig1 = sig*ones(1, r); % columns repeat
    sig1t = sig1'; % rows repeat
    s1 = sig1(:);
    s2 = sig1t(:);
    
    % The block elements
    a =  s1.^2 + s2.^2; % a column vector
    
    % Solve the linear system of equations
    cu = b1./a; %a.\b1;
    cv = b2./a; %a.\b2;
    
    % Matricize
    CU = reshape(cu, r, r);
    CV = reshape(cv, r, r);
    
    % Do the similarity transforms
    BU = U*CU*U';
    BV = V*CV*V';
    
    % %% Debug
    %
    % norm(R*R'*BU + BU*R*R' - E, 'fro');
    % norm((Sigma.^2)*CU + CU*(Sigma.^2) - E_mod, 'fro');
    % norm(a.*cu - b1, 'fro');
    %
    % norm(R'*R*BV + BV*R'*R - F, 'fro');
    %
    % BU1 = lyap(R*R', - E);
    % norm(R*R'*BU1 + BU1*R*R' - E, 'fro');
    %
    % BV1 = lyap(R'*R, - F);
    % norm(R'*R*BV1 + BV1*R'*R - F, 'fro');
    %
    % % as accurate as the lyap
    % norm(BU - BU1, 'fro')
    % norm(BV - BV1, 'fro')
end
        
    end