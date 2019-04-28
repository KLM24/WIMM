function L = my_bisection(Sigma, D, rho, bi_tol)
% this is right
%   Bisection algorithm - MATLAB implementation of Algorithm 1
%
%   Syntax: L = my_bisection(Sigma, D, rho, bi_tol)
%   my_bisection() computes the solution to the subproblem in the Frank-Wolfe algorithm.
%
%   Sigma:  Covariance matrix of the prior distribution
%   D:      The gradient of the objective function at an iteration
%   rho:    Wasserstein ambiguity size
%   bi_tol: Bisection tolerance value

d = size(D,1);

% Auxilary functions
h = @(inv_D) rho^2 - sum(sum(Sigma .* (eye(d) - inv_D)^2));
vec = @(x) x(:);

% Finding the bisection intervals
[v_1, lambda_1] = eigs(D,1);
LB = lambda_1 * (1 + sqrt(v_1'*Sigma*v_1)/rho);
UB = lambda_1 * (1 + sqrt(trace(Sigma))/rho);

% The main loop
while true
    gamma = (LB + UB)/2;
    
    D_inv = gamma * inv(gamma*eye(d) - D);
    L = D_inv * Sigma * D_inv;
    h_val = h(D_inv);
    
    if h_val < 0
        LB = gamma;
    else
        UB = gamma;
    end
    
    Delta = gamma * (rho^2 - trace(Sigma)) + gamma * vec(D_inv)'*vec(Sigma) - vec(L)'*vec(D);
    
    if (h_val >= 0) && (Delta < bi_tol)
        break
    end
end
end

