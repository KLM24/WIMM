function [phi_star, Q_star] = F_W(mu, Sigma, rho, x_dim)
% this is right
%   Frank-Wolfe algorithm - MATLAB implementation of Algorithm 2
%
%   Syntax: [phi_star, Q_star, obj, res] = Frank_Wolfe(mu, Sigma, rho, x_dim, opts)
%   Frank_Wolfe() computes the Nash equilibrium for the minimax problem.
%
%   Frank_Wolfe has the following input arguments:
%   mu:         Mean vector of the prior distribution
%   Sigma:      Covariance matrix of the prior distribution
%   rho:        Wasserstein ambiguity size
%   x_dim:      Dimension of the unobserved random variable x
%   opts:       Structure contains the parameters for Frank-Wolfe algorithm
%       opts.iter_max:  Maximum number of iteration
%                       By default:: 1000
%       opts.tol:       Relative duality gap as stopping criteria
%                       By default:: 1e-4
%       opts.verbose:   Set to true to turn on verbose output
%                       by default:: false
%   Frank_Wolfe returns the following outputs:
%   phi_star:   A structure contains parameters of the optimal affine decision rule
%       phi_star.G:     Slope of the optimal decision rule
%       phi_star.g:     Intercept of the optimal decision rule
%   Q_star:    A structure contains parameters of the least favorable prior
%       Q_star.mu:      Mean vector of the least favorable prior
%       Q_star.Sigma:   Covariance matrix of the least favorable prior
%   obj:       Stored objective value in any iterations of Frank-Wolfe algorithm
%   res:       Stored dual optimality gap in any iterations of Frank-Wolfe algorithm

% Parameters
n = x_dim;
sigma_s = eigs(Sigma,n);
sigma_script = sigma_s(n);
sigma_bar = (rho+sqrt(trace(Sigma)))^2;
C_bar = 2*sigma_bar^4/sigma_script^3;
iter_max = 1000;
tol = 0.01;

% Auxilary functions
G_ = @(S) S(1:n, n+1:end)/(S(n+1:end, n+1:end));
grad_f_ = @(G) [eye(n), -G]' * [eye(n), -G];

% Initialize
S = Sigma;

% Return Bayes estimator if rho is zero
if rho == 0
    G = G_(Sigma);
    phi_star.G = G;
    phi_star.g = mu(1:n) - G * mu(n+1:end);
    Q_star.Sigma = Sigma;
    Q_star.mu = mu;
    return
end


% The main loop
for iter = 0 : iter_max
    alpha = 2 / (2 + iter);
    % Printing the objective values
    G = G_(S);
    
    % Computing partial derivitive & solve the algebraic equation
    D = grad_f_(G);
    epsilon = alpha*tol*C_bar;
    L = my_bisection(Sigma, D, rho, epsilon);
%     re1 = vec(S)'* vec(D);
    
    S = S + alpha * (L-S);
%     re2 = vec(S)' * vec(grad_f_(G_(S)));
%     if abs(re2-re1) / re1 < tol
%         break
%     end
end

phi_star.G = G_(S);
phi_star.g = mu(1:n) - G_(S) * mu(n+1:end);
Q_star.Sigma = S;
Q_star.mu = mu;

end