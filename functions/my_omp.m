function [x_hat,index_set] = my_omp(A, y, K)
% Orthogonal Matching Pursuit (OMP) algorithm
%
% Inputs:
%   A : Measurement matrix (MxN)
%   y : Measurement vector (Mx1)
%   K : Sparsity level (number of non-zero elements in x)
%
% Output:
%   x_hat : Estimated sparse signal (Nx1)

    % Initialization
    [M, N] = size(A);
    residual = y;              % Initial residual
    index_set = [];            % Support set
    x_hat = zeros(N, 1);       % Initialize solution

    for k = 1:K
        %  Find the index with maximum correlation
        correlations = abs(A' * residual);
        [~, new_index] = max(correlations);
        index_set = unique([index_set, new_index]);

        %  Solve least squares problem on current support
        A_restricted = A(:, index_set);
        x_restricted = A_restricted \ y;

        % Update residual
        residual = y - A_restricted * x_restricted;

        %  early stopping if residual is small
        if norm(residual) < 1e-6
            break;
        end
    end

    % Step 4: Construct full x_hat
    x_hat(index_set) = x_restricted;
end
